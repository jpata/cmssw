#include "RecoTauTag/RecoTau/interface/RecoTauBinnedIsolationPlugin.h"
#include "RecoTauTag/RecoTau/interface/RecoTauIsolationMasking.h"

namespace reco::tau {

template<class Extractor>
class RecoTauDiscriminationBinnedIsolationImpl
: public RecoTauDiscriminationBinnedIsolation {
  public:
    RecoTauDiscriminationBinnedIsolationImpl(const edm::ParameterSet& pset)
      :RecoTauDiscriminationBinnedIsolation(pset),extractor_(pset) {}
  std::vector<reco::CandidatePtr> extractIsoObjects(
        const reco::PFTauRef& tau) const override {
      return extractor_(tau);
    }
  private:
    Extractor extractor_;
};

} // end namespace reco::tau


// Methods to get the right kind of canddiates
namespace {

class TrackExtractor {
  public:
    TrackExtractor(const edm::ParameterSet& pset){};
    std::vector<reco::CandidatePtr> operator()(const reco::PFTauRef& tau) const {
      return tau->isolationChargedHadrCands();
    }
};

class ECALExtractor {
  public:
    ECALExtractor(const edm::ParameterSet& pset){};
    std::vector<reco::CandidatePtr> operator()(const reco::PFTauRef& tau) const {
      return tau->isolationGammaCands();
    }
};

class MaskedECALExtractor {
  public:
    MaskedECALExtractor(const edm::ParameterSet& pset)
      :mask_(pset.getParameter<edm::ParameterSet>("mask")){};
    std::vector<reco::CandidatePtr> operator()(const reco::PFTauRef& tau) const {
      std::vector<reco::CandidatePtr> output;
      reco::tau::RecoTauIsolationMasking::IsoMaskResult
        result = mask_.mask(*tau);
      output.reserve(result.gammas.size());
      for(auto const& gamma : result.gammas) {
        output.push_back(gamma);
      }
      return output;
    }
  private:
    reco::tau::RecoTauIsolationMasking mask_;
};

class HCALExtractor {
  public:
    HCALExtractor(const edm::ParameterSet& pset){};
    std::vector<reco::CandidatePtr> operator()(const reco::PFTauRef& tau) const {
      return tau->isolationNeutrHadrCands();
    }
};

class MaskedHCALExtractor {
  public:
    MaskedHCALExtractor(const edm::ParameterSet& pset)
      :mask_(pset.getParameter<edm::ParameterSet>("mask")){};
    std::vector<reco::CandidatePtr> operator()(const reco::PFTauRef& tau) const {
      std::vector<reco::CandidatePtr> output;
      reco::tau::RecoTauIsolationMasking::IsoMaskResult
        result = mask_.mask(*tau);
      output.reserve(result.h0s.size());
      for(auto const& h0 : result.h0s) {
        output.push_back(h0);
      }
      return output;
    }
  private:
    reco::tau::RecoTauIsolationMasking mask_;
};
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_EDM_PLUGIN(RecoTauDiscriminantPluginFactory,
    reco::tau::RecoTauDiscriminationBinnedIsolationImpl<TrackExtractor>,
    "RecoTauDiscriminationBinnedTrackIsolation");
DEFINE_EDM_PLUGIN(RecoTauDiscriminantPluginFactory,
    reco::tau::RecoTauDiscriminationBinnedIsolationImpl<ECALExtractor>,
    "RecoTauDiscriminationBinnedECALIsolation");
DEFINE_EDM_PLUGIN(RecoTauDiscriminantPluginFactory,
    reco::tau::RecoTauDiscriminationBinnedIsolationImpl<MaskedECALExtractor>,
    "RecoTauDiscriminationBinnedMaskedECALIsolation");
DEFINE_EDM_PLUGIN(RecoTauDiscriminantPluginFactory,
    reco::tau::RecoTauDiscriminationBinnedIsolationImpl<HCALExtractor>,
    "RecoTauDiscriminationBinnedHCALIsolation");
DEFINE_EDM_PLUGIN(RecoTauDiscriminantPluginFactory,
    reco::tau::RecoTauDiscriminationBinnedIsolationImpl<MaskedHCALExtractor>,
    "RecoTauDiscriminationBinnedMaskedHCALIsolation");
