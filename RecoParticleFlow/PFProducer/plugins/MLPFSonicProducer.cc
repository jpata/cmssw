#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonEDProducer.h"
#include "RecoParticleFlow/PFProducer/interface/MLPFModel.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlockElementTrack.h"

//using namespace cms::Ort;

//use this to switch on detailed print statements in MLPF
//#define MLPF_DEBUG

class MLPFSonicProducer : public TritonEDProducer<> {
public:
  explicit MLPFSonicProducer(const edm::ParameterSet&);
  void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) override;
  void produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDPutTokenT<reco::PFCandidateCollection> pfCandidatesPutToken_;
  const edm::EDGetTokenT<reco::PFBlockCollection> inputTagBlocks_;
  unsigned int num_elements_total_;
  std::vector<const reco::PFBlockElement*> selected_elements_;
};

MLPFSonicProducer::MLPFSonicProducer(const edm::ParameterSet& cfg)
    : TritonEDProducer<>(cfg, "MLPFSonicProducer"),
      pfCandidatesPutToken_{produces<reco::PFCandidateCollection>()},
      inputTagBlocks_(consumes<reco::PFBlockCollection>(cfg.getParameter<edm::InputTag>("src"))),
      num_elements_total_(0) {}

void MLPFSonicProducer::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) {
  // one event per batch
  client_->setBatchSize(1);

  using namespace reco::mlpf;

  const auto& blocks = iEvent.get(inputTagBlocks_);
  const auto& all_elements = getPFElements(blocks);

  selected_elements_.clear();
  num_elements_total_ = 0;
  //unsigned int num_elements_total = 0;
  for (const auto* pelem : all_elements) {
    if (pelem->type() == reco::PFBlockElement::PS1 || pelem->type() == reco::PFBlockElement::PS2) {
      continue;
    }
    num_elements_total_ += 1;
    selected_elements_.push_back(pelem);
  }
  assert(num_elements_total_ < NUM_MAX_ELEMENTS_BATCH);

  //tensor size must be a multiple of the bin size and larger than the number of elements
  const auto tensor_size = LSH_BIN_SIZE * std::max(2u, (num_elements_total_ / LSH_BIN_SIZE + 1));
  assert(tensor_size <= NUM_MAX_ELEMENTS_BATCH);
  assert(tensor_size % LSH_BIN_SIZE == 0);

#ifdef MLPF_DEBUG
  std::cout << "tensor_size=" << tensor_size << std::endl;
#endif

  auto& input = iInput.at("x:0");
  input.setShape(0, tensor_size);
  auto inputdata = input.allocate<float>();
  auto& vinputdata = (*inputdata)[0];

  unsigned int ielem = 0;
  for (const auto* pelem : selected_elements_) {
    if (ielem > tensor_size) {
      //continue;
      break;
    }

    const auto& elem = *pelem;

    //prepare the input array from the PFElement
    const auto& props = getElementProperties(elem);

    //copy features to the input array
    for (unsigned int iprop = 0; iprop < NUM_ELEMENT_FEATURES; iprop++) {
      vinputdata.push_back(normalize(props[iprop]));
    }
    ielem += 1;
  }

  vinputdata.resize(NUM_ELEMENT_FEATURES * tensor_size);

  input.toServer(inputdata);
}

void MLPFSonicProducer::produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) {
  using namespace reco::mlpf;
  const auto& output1 = iOutput.begin()->second;
  const auto& outputs = output1.fromServer<float>();

  std::vector<reco::PFCandidate> pOutputCandidateCollection;
  for (size_t ielem = 0; ielem < num_elements_total_; ielem++) {
    std::vector<float> pred_id_probas(IDX_CLASS + 1, 0.0);
    const reco::PFBlockElement* elem = selected_elements_[ielem];

    for (unsigned int idx_id = 0; idx_id <= IDX_CLASS; idx_id++) {
      auto pred_proba = outputs[0][ielem * NUM_OUTPUT_FEATURES + idx_id];
      assert(!std::isnan(pred_proba));
      pred_id_probas[idx_id] = pred_proba;
    }

    auto imax = argMax(pred_id_probas);

    //get the most probable class PDGID
    int pred_pid = pdgid_encoding[imax];

    //a particle was predicted for this PFElement, otherwise it was a spectator
    if (pred_pid != 0) {
      //muons and charged hadrons should only come from tracks, otherwise we won't have track references to pass downstream
      if (((pred_pid == 13) || (pred_pid == 211)) && elem->type() != reco::PFBlockElement::TRACK) {
        pred_pid = 130;
      }

      if (elem->type() == reco::PFBlockElement::TRACK) {
        const auto* eltTrack = dynamic_cast<const reco::PFBlockElementTrack*>(elem);

        //a track with no muon ref should not produce a muon candidate, instead we interpret it as a charged hadron
        if (pred_pid == 13 && eltTrack->muonRef().isNull()) {
          pred_pid = 211;
        }

        //tracks from displaced vertices need reference debugging downstream as well, so we just treat them as neutrals for the moment
        if ((pred_pid == 211) && (eltTrack->isLinkedToDisplacedVertex())) {
          pred_pid = 130;
        }
      }

      //get the predicted momentum components
      float pred_pt = outputs[0][ielem * NUM_OUTPUT_FEATURES + IDX_PT];
      float pred_eta = outputs[0][ielem * NUM_OUTPUT_FEATURES + IDX_ETA];
      float pred_sin_phi = outputs[0][ielem * NUM_OUTPUT_FEATURES + IDX_SIN_PHI];
      float pred_cos_phi = outputs[0][ielem * NUM_OUTPUT_FEATURES + IDX_COS_PHI];
      float pred_e = outputs[0][ielem * NUM_OUTPUT_FEATURES + IDX_ENERGY];
      float pred_charge = outputs[0][ielem * NUM_OUTPUT_FEATURES + IDX_CHARGE];

      auto cand = makeCandidate(pred_pid, pred_charge, pred_pt, pred_eta, pred_sin_phi, pred_cos_phi, pred_e);
      setCandidateRefs(cand, selected_elements_, ielem);
      pOutputCandidateCollection.push_back(cand);

#ifdef MLPF_DEBUG
      std::cout << "ielem=" << ielem << " cand: pid=" << cand.pdgId() << " E=" << cand.energy() << " pt=" << cand.pt()
                << " eta=" << cand.eta() << " phi=" << cand.phi() << " charge=" << cand.charge() << std::endl;
#endif
    }
  }  //loop over PFElements

  iEvent.emplace(pfCandidatesPutToken_, pOutputCandidateCollection);
}

void MLPFSonicProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  TritonClient::fillPSetDescription(desc);
  desc.add<edm::InputTag>("src", edm::InputTag("particleFlowBlock"));
  //descriptions.addWithDefaultLabel(desc);
  descriptions.add("MLPFSonicProducer", desc);
  std::vector<const reco::PFBlockElement*> selected_elements;
}

DEFINE_FWK_MODULE(MLPFSonicProducer);
