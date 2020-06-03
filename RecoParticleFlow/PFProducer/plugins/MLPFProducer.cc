#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementGsfTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementBrem.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#include "TLorentzVector.h"
#include "TStopwatch.h"

//The model takes the following number of features for each input PFElement
static const unsigned int NUM_ELEMENT_FEATURES = 15;

//The model has 13 outputs for each particle, which consist of:

static const unsigned int NUM_OUTPUTS = 12;
//0...7: class probabilities
//8: eta
//9: phi
//10: energy
//11: charge
static const unsigned int NUM_CLASS = 7;
static const unsigned int IDX_ETA = 8;
static const unsigned int IDX_PHI = 9;
static const unsigned int IDX_ENERGY = 10;
static const unsigned int IDX_CHARGE = 11;


//index [0, N_pdgids) -> PDGID
//this maps the absolute values of the predicted PDGIDs to an array of ascending indices
static const std::vector<int> pdgid_encoding = {0, 1, 2, 11, 13, 22, 130, 211};

//PFElement::type -> index [0, N_types)
//this maps the type of the PFElement to an ascending index that is used by the model to distinguish between different elements
static const std::map<int, int> elem_type_encoding = {
  {1, 0},
  {2, 1},
  {3, 2},
  {4, 3},
  {5, 4},
  {6, 5},
  {7, 6},
  {8, 7},
  {9, 8},
  {10, 9},
  {11, 10},
};

struct MLPFCache {
  std::atomic<tensorflow::GraphDef*> graph_def;
};

//Prepares the input array of floats for a single PFElement 
std::array<float, NUM_ELEMENT_FEATURES> get_element_properties(const reco::PFBlockElement& orig) {
  const auto type = orig.type();
  float pt = 0.0;
  //these are placeholders for the the future
  [[maybe_unused]] float deltap = 0.0;
  [[maybe_unused]] float sigmadeltap = 0.0;
  [[maybe_unused]] float px = 0.0;
  [[maybe_unused]] float py = 0.0;
  [[maybe_unused]] float pz = 0.0;
  float eta = 0.0;
  float phi = 0.0;
  float energy = 0.0;
  float trajpoint = 0.0;
  float eta_ecal = 0.0;
  float phi_ecal = 0.0;
  float eta_hcal = 0.0;
  float phi_hcal = 0.0;
  float charge = 0;
  float layer = 0;
  float depth = 0;
  float muon_dt_hits  = 0.0;
  float muon_csc_hits  = 0.0;

  if (type == reco::PFBlockElement::TRACK) {
    const auto& matched_pftrack = orig.trackRefPF();
    if (matched_pftrack.isNonnull()) {
      const auto& atECAL = matched_pftrack->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
      const auto& atHCAL = matched_pftrack->extrapolatedPoint(reco::PFTrajectoryPoint::HCALEntrance);
      if (atHCAL.isValid()) {
        eta_hcal = atHCAL.positionREP().eta();
        phi_hcal = atHCAL.positionREP().phi();
      }
      if (atECAL.isValid()) {
        eta_ecal = atECAL.positionREP().eta();
        phi_ecal = atECAL.positionREP().phi();
      }
    }
    const auto& ref = ((const reco::PFBlockElementTrack*)&orig)->trackRef();
    pt = ref->pt();
    px = ref->px();
    py = ref->py();
    pz = ref->pz();
    eta = ref->eta();
    phi = ref->phi();
    energy = ref->pt() * cosh(ref->eta());
    charge = ref->charge();

    reco::MuonRef muonRef = orig.muonRef();
    if (muonRef.isNonnull()) {
      reco::TrackRef standAloneMu = muonRef->standAloneMuon();
      if (standAloneMu.isNonnull()) {
        muon_dt_hits = standAloneMu->hitPattern().numberOfValidMuonDTHits();
        muon_csc_hits = standAloneMu->hitPattern().numberOfValidMuonCSCHits();
      }
    }

  } else if (type == reco::PFBlockElement::BREM) {
    const auto* orig2 = (const reco::PFBlockElementBrem*)&orig;
    const auto& ref = orig2->GsftrackRef();
    if (ref.isNonnull()) {
      deltap = orig2->DeltaP();
      sigmadeltap = orig2->SigmaDeltaP();
      pt = ref->pt();
      px = ref->px();
      py = ref->py();
      pz = ref->pz();
      eta = ref->eta();
      phi = ref->phi();
      energy = ref->pt() * cosh(ref->eta());
      trajpoint = orig2->indTrajPoint();
      charge = ref->charge();
    }
  } else if (type == reco::PFBlockElement::GSF) {
    //requires to keep GsfPFRecTracks
    const auto* orig2 = (const reco::PFBlockElementGsfTrack*)&orig;
    pt = orig2->Pin().pt();
    px = orig2->Pin().px();
    py = orig2->Pin().py();
    pz = orig2->Pin().pz();
    eta = orig2->Pin().eta();
    phi = orig2->Pin().phi();
    energy = pt * cosh(eta);
    if (!orig2->GsftrackRefPF().isNull()) {
      charge = orig2->GsftrackRefPF()->charge();
    }
  } else if (type == reco::PFBlockElement::ECAL || type == reco::PFBlockElement::PS1 ||
             type == reco::PFBlockElement::PS2 || type == reco::PFBlockElement::HCAL ||
             type == reco::PFBlockElement::HO ||
             type == reco::PFBlockElement::HFHAD || type == reco::PFBlockElement::HFEM) {
    const auto& ref = ((const reco::PFBlockElementCluster*)&orig)->clusterRef();
    if (ref.isNonnull()) {
      eta = ref->eta();
      phi = ref->phi();
      px = ref->position().x();
      py = ref->position().y();
      pz = ref->position().z();
      energy = ref->energy();
      layer = ref->layer();
      depth = ref->depth();
    }
  } else if (type == reco::PFBlockElement::SC) {
    const auto& clref = ((const reco::PFBlockElementSuperCluster*)&orig)->superClusterRef();
    if (clref.isNonnull()) {
      eta = clref->eta();
      phi = clref->phi();
      px = clref->position().x();
      py = clref->position().y();
      pz = clref->position().z();
      energy = clref->energy();
    }
  }

  float typ_idx = static_cast<float>(elem_type_encoding.at(orig.type()));
 
  //Must be the same order as in tf_model.py 
  return std::array<float, NUM_ELEMENT_FEATURES>({{
    typ_idx,
    pt, eta, phi, energy,
    layer, depth, charge, trajpoint,
    eta_ecal, phi_ecal, eta_hcal, phi_hcal,
    muon_dt_hits, muon_csc_hits}});

}

class MLPFProducer : public edm::stream::EDProducer<edm::GlobalCache<MLPFCache> > {
public:
  explicit MLPFProducer(const edm::ParameterSet&, const MLPFCache*);
  void produce(edm::Event& event, const edm::EventSetup& setup) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  // static methods for handling the global cache
  static std::unique_ptr<MLPFCache> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(MLPFCache*);

private:
  const edm::EDPutTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
  const edm::EDGetTokenT<reco::PFBlockCollection> inputTagBlocks_;
  const unsigned int min_batch_size_;
  const std::string model_path_;
  tensorflow::Session* session_;
  TStopwatch sw;
};

MLPFProducer::MLPFProducer(const edm::ParameterSet& cfg, const MLPFCache* cache) :
  pfCandidatesToken_{produces<reco::PFCandidateCollection>()},
  inputTagBlocks_(consumes<reco::PFBlockCollection>(cfg.getParameter<edm::InputTag>("src"))),
  min_batch_size_(cfg.getParameter<unsigned int>("min_batch_size")),
  model_path_(cfg.getParameter<std::string>("model_path")) {
    session_ = tensorflow::createSession(cache->graph_def);
}

//returns the argmax of a given std::vector<T>
template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

float normalize(float in) {
    if (std::abs(in) > 1e4) {
        return 0.0;
    }
    else if (std::isnan(in)) {
        return 0.0;
    }
    return in;
}

void MLPFProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  sw.Reset();

  auto blocks_handle = event.getHandle(inputTagBlocks_);
  const auto& blocks = *blocks_handle;
  std::vector<reco::PFCandidate> pOutputCandidateCollection;

  //Create input batches from PFElements in the event.
  //Elements can only pass information to each other within the same batch,
  //similar to how PFAlgo processes each PFBlock independently.
  std::vector<std::vector<const reco::PFBlockElement*>> element_batches;
  std::vector<const reco::PFBlockElement*> element_batch;
  for (const auto& block : blocks) {
    const auto& elems = block.elements();
    if (element_batch.size () > min_batch_size_) {
      element_batches.push_back(element_batch);
      element_batch.clear();
    }
    for (const auto& elem : elems) {
      element_batch.push_back(&elem);
    }
  }
  //last chunk
  element_batches.push_back(element_batch);
  element_batch.clear();

  //get the maximum number of elements for all the batches
  unsigned int num_max_elems_batch = 0;
  unsigned int num_elements_total = 0;
 
  std::cout << "num_batches=" << element_batches.size();
  for (const auto& element_batch : element_batches) {
    if (element_batch.size() > num_max_elems_batch) {
      num_max_elems_batch = element_batch.size();
    }
    std::cout << " " << element_batch.size();
    num_elements_total += element_batch.size();
  }
  std::cout << " num_max_elems_batch=" << num_max_elems_batch << std::endl;

  const unsigned int num_batches = element_batches.size();

  //Create the input tensor
  tensorflow::TensorShape shape({num_batches, num_max_elems_batch, NUM_ELEMENT_FEATURES});
  tensorflow::Tensor input(tensorflow::DT_FLOAT, shape);
  input.flat<float>().setZero();
    
  //Fill the input tensor
  sw.Start();
  unsigned int ibatch = 0;
  for (const auto& element_batch : element_batches) {

    //Put all PFElements in this PFBlock into the input tensor.
    unsigned int ielem = 0;
    for (const auto* pelem : element_batch) {

      //get the PFElement as a pointer
      const auto& elem = *pelem;

      //prepare the input array from the PFElement
      const auto& props = get_element_properties(elem);

      //copy features to the input array
      for (unsigned int iprop=0; iprop<NUM_ELEMENT_FEATURES; iprop++) {
        input.tensor<float, 3>()(ibatch, ielem, iprop) = normalize(props[iprop]);
      }
      ielem += 1; 
    }

    ibatch += 1;
  } //loop over PFBlocks
  tensorflow::NamedTensorList input_list = {{"x:0", input}};

  //Prepare the output tensor
  std::vector<tensorflow::Tensor> outputs;
  std::vector<std::string> output_names = {"Identity:0"};

  //run the GNN inference, given the inputs and the output.
  //Note that the GNN enables information transfer between the input PFElements,
  //such that the output ML-PFCandidates are in general combinations of the input PFElements, in the form of
  //y_out = Adj.x_in, where x_in is input matrix (num_elem, NUM_ELEMENT_FEATURES), y_out is the output matrix (num_elem, NUM_OUTPUT_FEATURES)
  //and Adj is an adjacency matrix between the elements that is constructed on the fly during model inference.
  tensorflow::run(session_, input_list, output_names, &outputs);

  //process the output tensor to ML-PFCandidates.
  //The output can contain up to num_elem particles, with predicted PDGID=0 corresponding to no particles predicted.
  const auto out_arr = outputs[0].tensor<float, 3>();

  unsigned int num_pred_particles = 0;
  for (unsigned int ibatch=0; ibatch<num_batches; ibatch++) {
    for (unsigned int ielem=0; ielem<element_batches.at(ibatch).size(); ielem++) {

      //get the coefficients in the output corresponding to the class probabilities (raw logits) 
      std::vector<float> pred_id_logits;
      for (unsigned int idx_id=0; idx_id <= NUM_CLASS; idx_id++) {
        pred_id_logits.push_back(out_arr(ibatch, ielem, idx_id));
      }
      //get the most probable class PDGID
      int pred_pid = pdgid_encoding.at(arg_max(pred_id_logits));

      //get the predicted momentum components
      float pred_eta = out_arr(ibatch, ielem, IDX_ETA);
      float pred_phi = out_arr(ibatch, ielem, IDX_PHI);
      pred_phi = TMath::ATan2(TMath::Sin(pred_phi), TMath::Cos(pred_phi));
      float pred_e = out_arr(ibatch, ielem, IDX_ENERGY);

      //currently, set the pT from a massless approximation.
      //later versions of the model may predict predict both the energy and pT of the particle
      float pred_pt = pred_e / cosh(pred_eta);

      float pred_charge = out_arr(ibatch, ielem, IDX_CHARGE);

      //very verbose
      //element_batches.at(ibatch).at(ielem)->Dump();
      std::cout << "MLPF " << pred_pid << "," << pred_pt << "," << pred_e << "," << pred_eta << "," << pred_phi << std::endl;

      //a particle was predicted for this PFElement
      if (pred_pid != 0) {

        //set the charge to +1 or -1 for PFCandidates that are charged, according to the sign of the predicted charge
        reco::PFCandidate::Charge charge = 0;
        if (pred_pid == 11 || pred_pid == 13 || pred_pid == 211) {
            charge = pred_charge > 0 ? +1 : -1;
        }

        TLorentzVector p4;
        p4.SetPtEtaPhiE(pred_pt, pred_eta, pred_phi, pred_e);

        reco::PFCandidate cand(0, math::XYZTLorentzVector(p4.X(), p4.Y(), p4.Z(), p4.E()), reco::PFCandidate::ParticleType(0));
        cand.setPdgId(pred_pid);
        cand.setCharge(charge);

        const reco::PFBlockElement* elem = element_batches.at(ibatch).at(ielem);
        if (elem->type() == reco::PFBlockElement::TRACK) {
            cand.setTrackRef(elem->trackRef());
        }

        pOutputCandidateCollection.push_back(cand);
        num_pred_particles += 1;
      } else {
        //this element did not directly yield a particle, but could have been used indirectly for other ML-PFCandidates
      }
    } //loop over PFElements
  }
  sw.Stop();

  std::cout << "num_elements_total=" << num_elements_total
    << " num_pred_particles=" << num_pred_particles
    << " cputime=" << sw.CpuTime() << std::endl;

  event.emplace(pfCandidatesToken_, pOutputCandidateCollection);
}

std::unique_ptr<MLPFCache> MLPFProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  // this method is supposed to create, initialize and return a MLPFCache instance
  std::unique_ptr<MLPFCache> cache = std::make_unique<MLPFCache>();

  //load the frozen TF graph of the GNN model
  std::string path = params.getParameter<std::string>("model_path");
  auto fullPath = edm::FileInPath(path).fullPath();
  cache->graph_def = tensorflow::loadGraphDef(fullPath);

  return cache;
}

void MLPFProducer::globalEndJob(MLPFCache* cache) {
  delete cache->graph_def;
}

void MLPFProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  desc.add<unsigned int>("min_batch_size", 500);
  desc.add<std::string>("model_path", "RecoParticleFlow/PFProducer/data/mlpf/mlpf_2020_05_07.pb");
  descriptions.add("MLPFProducer", desc);
}

DEFINE_FWK_MODULE(MLPFProducer);
