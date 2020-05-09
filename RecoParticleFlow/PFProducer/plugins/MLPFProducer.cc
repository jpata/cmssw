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

static const unsigned int NUM_ELEMENT_FEATURES = 15;

//index [0, N_pdgids) -> PDGID
static const std::vector<int> pdgid_encoding = {0, 1, 2, 11, 13, 22, 130, 211};

//PFElement::type -> index [0, N_types)
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

std::array<float, NUM_ELEMENT_FEATURES> get_element_properties(const reco::PFBlockElement& orig) {
  const auto type = orig.type();
  float pt = 0.0;
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
    typ_idx, pt, eta, phi,
    energy, layer, depth, charge,
    trajpoint, eta_ecal, phi_ecal, eta_hcal, phi_hcal,
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
  const std::string model_path_;
  tensorflow::Session* session_;
};

MLPFProducer::MLPFProducer(const edm::ParameterSet& cfg, const MLPFCache* cache) :
  pfCandidatesToken_{produces<reco::PFCandidateCollection>()},
  inputTagBlocks_(consumes<reco::PFBlockCollection>(cfg.getParameter<edm::InputTag>("src"))),
  model_path_(cfg.getParameter<std::string>("model_path")) {
    session_ = tensorflow::createSession(cache->graph_def);
}


template <typename T, typename A>
int arg_max(std::vector<T, A> const& vec) {
  return static_cast<int>(std::distance(vec.begin(), max_element(vec.begin(), vec.end())));
}

void MLPFProducer::produce(edm::Event& event, const edm::EventSetup& setup) {
  auto blocks_handle = event.getHandle(inputTagBlocks_);
  const auto& blocks = *blocks_handle;
  std::vector<reco::PFCandidate> pOutputCandidateCollection;

  //Put all PFElements into the input tensor
  for (const auto& block : blocks) {
    const auto& elems = block.elements();
    unsigned int num_elems = elems.size();
    if (num_elems > 1) {
      std::cout << "num_elems=" << num_elems << std::endl;
    }

    tensorflow::TensorShape shape({num_elems, NUM_ELEMENT_FEATURES});
    tensorflow::Tensor input(tensorflow::DT_FLOAT, shape);
    input.flat<float>().setZero();

    //create input array from PFElements
    unsigned int ielem = 0;
    for (const auto& elem : elems) {
      const auto& props = get_element_properties(elem);
      for (unsigned int iprop=0; iprop<NUM_ELEMENT_FEATURES; iprop++) {
        input.tensor<float, 2>()(ielem, iprop) = props[iprop];
      }
      ielem += 1; 
    }

    tensorflow::NamedTensorList input_list = {{"x:0", input}};
    std::vector<tensorflow::Tensor> outputs;
    std::vector<std::string> output_names = {"Identity:0"};

    //run the GNN inference
    tensorflow::run(session_, input_list, output_names, &outputs);

    //process the outputs
    const auto out_arr = outputs[0].tensor<float, 2>();
    for (unsigned int ielem=0; ielem<num_elems; ielem++) {

      //get the coefficients in the output corresponding to the class probabilities (raw logits) 
      std::vector<float> pred_id_logits;
      for (unsigned int idx_id=0; idx_id<8; idx_id++) {
        pred_id_logits.push_back(out_arr(ielem, idx_id));
      }
      //get the most probable class PDGID
      int pred_pid = pdgid_encoding.at(arg_max(pred_id_logits));

      //get the predicted momentum components
      float pred_eta = out_arr(ielem, 8);
      float pred_phi = out_arr(ielem, 9);
      float pred_e = out_arr(ielem, 10);

      //massless approximation
      float pred_pt = pred_e / TMath::CosH(pred_eta);

      float pred_charge = out_arr(ielem, 11);

      std::cout << elems[ielem] << std::endl;
      //a particle was predicted for this PFElement
      if (pred_pid != 0) {

        reco::PFCandidate::Charge charge = 0;
        if (pred_pid == 11 || pred_pid == 13 || pred_pid == 211) {
            charge = pred_charge > 0 ? +1 : -1;
        }

        std::cout << "MLPFCandidate " <<
           " pid=" << pred_pid << " pt=" << pred_pt << " eta=" << pred_eta << " phi=" << pred_phi << " e=" << pred_e <<
           " charge=" << charge << std::endl;

        TLorentzVector p4;
        p4.SetPtEtaPhiE(pred_pt, pred_eta, pred_phi, pred_e);

        reco::PFCandidate cand(0, math::XYZTLorentzVector(p4.X(), p4.Y(), p4.Z(), p4.E()), reco::PFCandidate::ParticleType(0));
        cand.setParticleType(cand.translatePdgIdToType(pred_pid));
        cand.setCharge(charge);
        pOutputCandidateCollection.push_back(cand);
      } else {
        //this element did not directly yield a particle, but could have been used indirectly for other ML-PFCandidates
      }
    }

  }

  event.emplace(pfCandidatesToken_, pOutputCandidateCollection);
}

std::unique_ptr<MLPFCache> MLPFProducer::initializeGlobalCache(const edm::ParameterSet& params) {
  // this method is supposed to create, initialize and return a MLPFCache instance
  std::unique_ptr<MLPFCache> cache = std::make_unique<MLPFCache>();

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
  //desc.add<bool>("ignore_leptons", false);
  //desc.add<double>("norm_factor", 50.);
  //desc.add<unsigned int>("max_n_pf", 4500);
  desc.add<std::string>("model_path", "RecoParticleFlow/PFProducer/data/mlpf/mlpf_2020_05_07.pb");
  descriptions.add("MLPFProducer", desc);
}

DEFINE_FWK_MODULE(MLPFProducer);
