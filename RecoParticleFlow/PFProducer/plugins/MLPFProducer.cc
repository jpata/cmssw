#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "RecoParticleFlow/PFProducer/interface/MLPFModel.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"

struct MLPFCache {
  const tensorflow::GraphDef* graph_def;
};

class MLPFProducer : public edm::stream::EDProducer<edm::GlobalCache<MLPFCache> > {
public:
  explicit MLPFProducer(const edm::ParameterSet&, const MLPFCache*);
  void produce(edm::Event& event, const edm::EventSetup& setup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // static methods for handling the global cache
  static std::unique_ptr<MLPFCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(MLPFCache*);

private:
  const edm::EDPutTokenT<reco::PFCandidateCollection> pfCandidatesPutToken_;
  const edm::EDGetTokenT<reco::PFBlockCollection> inputTagBlocks_;
  const std::string model_path_;
  tensorflow::Session* session_;
};

MLPFProducer::MLPFProducer(const edm::ParameterSet& cfg, const MLPFCache* cache)
    : pfCandidatesPutToken_{produces<reco::PFCandidateCollection>()},
      inputTagBlocks_(consumes<reco::PFBlockCollection>(cfg.getParameter<edm::InputTag>("src"))),
      model_path_(cfg.getParameter<std::string>("model_path")) {
  session_ = tensorflow::createSession(cache->graph_def);
}

void MLPFProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  using namespace reco::mlpf;

  const auto& blocks = event.get(inputTagBlocks_);
  const auto& all_elements = getPFElements(blocks);

  const long long int num_elements_total = all_elements.size();

  //tensor size must be a multiple of the bin size and larger than the number of elements
  const auto tensor_size = LSH_BIN_SIZE * (num_elements_total / LSH_BIN_SIZE + 1);
  //const auto tensor_size = 6400;
  assert(tensor_size <= NUM_MAX_ELEMENTS_BATCH);
  assert(tensor_size % LSH_BIN_SIZE == 0);
  //std::cout << "tensor_size=" << tensor_size << std::endl;

  //Create the input tensor
  tensorflow::TensorShape shape({BATCH_SIZE, tensor_size, NUM_ELEMENT_FEATURES});
  tensorflow::Tensor input(tensorflow::DT_FLOAT, shape);
  input.flat<float>().setZero();

  //Fill the input tensor
  unsigned int ielem = 0;
  for (const auto* pelem : all_elements) {
    if (ielem > tensor_size) {
      continue;
    }

    const auto& elem = *pelem;

    // if (elem.type() == reco::PFBlockElement::PS1 || elem.type() == reco::PFBlockElement::PS2) {
    //   continue;
    // }

    //prepare the input array from the PFElement
    const auto& props = getElementProperties(elem);

    //copy features to the input array
    for (unsigned int iprop = 0; iprop < NUM_ELEMENT_FEATURES; iprop++) {
      input.tensor<float, 3>()(0, ielem, iprop) = normalize(props[iprop]);
    }
    ielem += 1;
  }
  //std::cout << "ielem=" << ielem << std::endl;

  //TF model input and output tensor names
  const tensorflow::NamedTensorList input_list = {{"x:0", input}};
  const std::vector<std::string> output_names = {"Identity:0"};

  //Prepare the output tensor
  std::vector<tensorflow::Tensor> outputs;

  //run the GNN inference, given the inputs and the output.
  //Note that the GNN enables information transfer between the input PFElements,
  //such that the output ML-PFCandidates are in general combinations of the input PFElements, in the form of
  //y_out = Adj.x_in, where x_in is input matrix (num_elem, NUM_ELEMENT_FEATURES), y_out is the output matrix (num_elem, NUM_OUTPUT_FEATURES)
  //and Adj is an adjacency matrix between the elements that is constructed on the fly during model inference.
  tensorflow::run(session_, input_list, output_names, &outputs);

  //process the output tensor to ML-PFCandidates.
  //The output can contain up to num_elem particles, with predicted PDGID=0 corresponding to no particles predicted.
  const auto out_arr = outputs[0].tensor<float, 3>();

  std::vector<reco::PFCandidate> pOutputCandidateCollection;
  for (unsigned int ielem = 0; ielem < all_elements.size(); ielem++) {
    //get the coefficients in the output corresponding to the class probabilities (raw logits)
    std::vector<float> pred_id_logits;
    for (unsigned int idx_id = 0; idx_id <= IDX_CLASS; idx_id++) {
      auto pred_logit = out_arr(0, ielem, idx_id);
      assert(!std::isnan(pred_logit));
      pred_id_logits.push_back(pred_logit);
      //std::cout << pred_logit << " ";
    }

    //get the most probable class PDGID
    int pred_pid = pdgid_encoding[argMax(pred_id_logits)];

    // if the charged candidate is generated from a track linked to a displaced vertex, we get
    // JetTracksAssociatorExplicit/'ak4JetTracksAssociatorExplicitAll' -> Ref is inconsistent with RefVectorid
    const reco::PFBlockElement* elem = all_elements[ielem];
    //std::cout << " pid=" << pred_pid << " type=" << elem->type() << std::endl;

    //at the moment, I don't know what kind of specific configuration of refs & other metadata the muons expect downstream of PF
    //so I currently reconstruct them as charged hadrons as a workaround to avoid crashes downstream.
    //this needs debugging of the downstream reco modules to understand what defined a PFCandidate muon, besides the 4-momentum
    if (pred_pid == 13) {
      pred_pid = 211;
    }

    //muons and charged hadrons should only come from tracks, otherwise we won't have references to pass downstream
    if (((pred_pid == 13) || (pred_pid==211)) && elem->type()!=reco::PFBlockElement::TRACK) {
      pred_pid = 130;
    }

    if (elem->type()==reco::PFBlockElement::TRACK) {
      const auto* eltTrack = dynamic_cast<const reco::PFBlockElementTrack*>(elem);

      //a track with no muon ref should not produce a muon candidate, instead we interpret it as a charged hadron
      if (pred_pid == 13 && eltTrack->muonRef().isNull()) {
        pred_pid = 211;
      }

      //tracks from displaced vertices need reference debugging downstream as well, so we just treat them as neutrals for the moment
      if ((pred_pid==211) && (eltTrack->isLinkedToDisplacedVertex())) {
        pred_pid = 130;
      }
    }

    //get the predicted momentum components
    float pred_pt = out_arr(0, ielem, IDX_PT);
    float pred_eta = out_arr(0, ielem, IDX_ETA);
    float pred_sin_phi = out_arr(0, ielem, IDX_SIN_PHI);
    float pred_cos_phi = out_arr(0, ielem, IDX_COS_PHI);
    float pred_e = out_arr(0, ielem, IDX_ENERGY);
    float pred_charge = out_arr(0, ielem, IDX_CHARGE);

    //a particle was predicted for this PFElement, otherwise it was a spectator
    if (pred_pid != 0) {
      //std::cout << pred_charge << " " << pred_pt << " " << pred_eta << " " << pred_sin_phi << " " << pred_cos_phi << " " << pred_e << std::endl;
      auto cand = makeCandidate(pred_pid, pred_charge, pred_pt, pred_eta, pred_sin_phi, pred_cos_phi, pred_e);
      setCandidateRefs(cand, all_elements, ielem);
      pOutputCandidateCollection.push_back(cand);
    }
  }  //loop over PFElements

  event.emplace(pfCandidatesPutToken_, pOutputCandidateCollection);
}

std::unique_ptr<MLPFCache> MLPFProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  // this method is supposed to create, initialize and return a MLPFCache instance
  std::unique_ptr<MLPFCache> cache = std::make_unique<MLPFCache>();

  //load the frozen TF graph of the GNN model
  std::string path = params.getParameter<std::string>("model_path");
  auto fullPath = edm::FileInPath(path).fullPath();
  LogDebug("MLPFProducer") << "Initializing MLPF model from " << fullPath;

  cache->graph_def = tensorflow::loadGraphDef(fullPath);

  return cache;
}

void MLPFProducer::globalEndJob(MLPFCache* cache) { delete cache->graph_def; }

void MLPFProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  desc.add<std::string>("model_path", "RecoParticleFlow/PFProducer/data/mlpf/mlpf_2021_02_12.pb");
  descriptions.addWithDefaultLabel(desc);
}

DEFINE_FWK_MODULE(MLPFProducer);
