#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerBase.h"
#include "CommonTools/RecoAlgos/interface/KDTreeLinkerAlgo.h"
#include "FWCore/SOA/interface/Column.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoParticleFlow/PFProducer/interface/Tables.h"

#include "TMath.h"

using namespace edm::soa;
using namespace edm::soa::col;

// This class is used to find all links between Tracks and ECAL clusters
// using a KDTree algorithm.
// It is used in PFBlockAlgo.cc in the function links().
class KDTreeLinkerTrackEcal : public KDTreeLinkerBase {
public:
  KDTreeLinkerTrackEcal(const edm::ParameterSet& conf);
  ~KDTreeLinkerTrackEcal() override;

  // With this method, we create the list of psCluster that we want to link.
  void insertTargetElt(reco::PFBlockElement* track) override;

  // Here, we create the list of ecalCluster that we want to link. From ecalCluster
  // and fraction, we will create a second list of rechits that will be used to
  // build the KDTree.
  void insertFieldClusterElt(reco::PFBlockElement* ecalCluster) override;

  // The KDTree building from rechits list.
  void buildTree(const PFTables& pftables) override;

  // Here we will iterate over all tracks. For each track intersection point with ECAL,
  // we will search the closest rechits in the KDTree, from rechits we will find the
  // ecalClusters and after that we will check the links between the track and
  // all closest ecalClusters.
  void searchLinks(const PFTables& pftables) override;

  // Here, we will store all PS/ECAL founded links in the PFBlockElement class
  // of each psCluster in the PFmultilinks field.
  void updatePFBlockEltWithLinks(const PFTables& pftables) override;

  // Here we free all allocated structures.
  void clear() override;

private:
  // Data used by the KDTree algorithm : sets of Tracks and ECAL clusters.
  BlockEltSet targetSet_;
  BlockEltSet fieldClusterSet_;
  std::vector<reco::PFBlockElement*> fieldClusterVec_;

  // Sets of rechits that compose the ECAL clusters.
  RecHitSet rechitsSet_;
  std::vector<const reco::PFRecHit*> rechitsVec_;

  // Map of linked Track/ECAL clusters.
  BlockElt2BlockEltMap target2ClusterLinks_;

  // Map of the ECAL clusters associated to a rechit.
  RecHit2BlockEltMap rechit2ClusterLinks_;
  std::unordered_map<size_t, std::vector<size_t>> rechit2ClusterLinksIdx_;

  // KD trees
  KDTreeLinkerAlgo<size_t> tree_;

  TrackTable trackTable_;
  ClusterTable clusterTable_;
  RecHitTable rechitTable_;
};

// the text name is different so that we can easily
// construct it when calling the factory
DEFINE_EDM_PLUGIN(KDTreeLinkerFactory, KDTreeLinkerTrackEcal, "KDTreeTrackAndECALLinker");

KDTreeLinkerTrackEcal::KDTreeLinkerTrackEcal(const edm::ParameterSet& conf) : KDTreeLinkerBase(conf) {}

KDTreeLinkerTrackEcal::~KDTreeLinkerTrackEcal() { clear(); }

void KDTreeLinkerTrackEcal::insertTargetElt(reco::PFBlockElement* track) {
  if (track->trackRefPF()->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax).isValid()) {
    targetSet_.insert(track);
  }
}

void KDTreeLinkerTrackEcal::insertFieldClusterElt(reco::PFBlockElement* ecalCluster) {
  const reco::PFClusterRef& clusterref = ecalCluster->clusterRef();

  // This test is more or less done in PFBlockAlgo.h. In others cases, it should be switch on.
  //   if (!((clusterref->layer() == PFLayer::ECAL_ENDCAP) ||
  // 	(clusterref->layer() == PFLayer::ECAL_BARREL)))
  //     return;

  const std::vector<reco::PFRecHitFraction>& fraction = clusterref->recHitFractions();

  // We create a list of ecalCluster
  fieldClusterSet_.insert(ecalCluster);
  for (size_t rhit = 0; rhit < fraction.size(); ++rhit) {
    const reco::PFRecHitRef& rh = fraction[rhit].recHitRef();
    double fract = fraction[rhit].fraction();

    if ((rh.isNull()) || (fract < cutOffFrac))
      continue;

    const reco::PFRecHit& rechit = *rh;

    // We save the links rechit to EcalClusters
    rechit2ClusterLinks_[&rechit].insert(ecalCluster);

    // We create a liste of rechits
    rechitsSet_.insert(&rechit);
  }
}

void KDTreeLinkerTrackEcal::buildTree(const PFTables& pftables) {
  //convert sets to ordered vectors
  for (const auto* rh : rechitsSet_) {
    rechitsVec_.push_back(rh);
  }
  rechitsSet_.clear();

  for (auto* cl : fieldClusterSet_) {
    fieldClusterVec_.push_back(cl);
  }
  fieldClusterSet_.clear();

  //trackTable_ = makeTrackTable(targetSet_);
  trackTable_ = pftables.track_table_;
  rechitTable_ = makeRecHitTable(rechitsVec_);
  clusterTable_ = makeClusterTable(fieldClusterVec_);

  //convert pointer-based links to index-based links
  for (const auto& rh_cluster : rechit2ClusterLinks_) {
    const auto* rechit = rh_cluster.first;
    const auto idx_rechit =
        std::distance(rechitsVec_.begin(), std::find(rechitsVec_.begin(), rechitsVec_.end(), rechit));
    for (const auto* cluster : rh_cluster.second) {
      const auto idx_cluster =
          std::distance(fieldClusterVec_.begin(), std::find(fieldClusterVec_.begin(), fieldClusterVec_.end(), cluster));
      rechit2ClusterLinksIdx_[idx_rechit].push_back(idx_cluster);
    }
  }
  rechit2ClusterLinks_.clear();

  // List of pseudo-rechits that will be used to create the KDTree
  std::vector<KDTreeNodeInfo<size_t, 2>> eltList;

  // Filling of this list
  for (size_t irechit = 0; irechit < rechitsVec_.size(); irechit++) {
    KDTreeNodeInfo<size_t, 2> rh1(
        irechit, float(rechitTable_.get<PF::rechit::Eta>(irechit)), float(rechitTable_.get<PF::rechit::Phi>(irechit)));
    eltList.push_back(rh1);

    // Here we solve the problem of phi circular set by duplicating some rechits
    // too close to -Pi (or to Pi) and adding (substracting) to them 2 * Pi.
    if (rh1.dims[1] > (M_PI - phiOffset_)) {
      float phi = rh1.dims[1] - 2 * M_PI;
      KDTreeNodeInfo<size_t, 2> rh2(irechit, float(rechitTable_.get<PF::rechit::Eta>(irechit)), phi);
      eltList.push_back(rh2);
    }

    if (rh1.dims[1] < (M_PI * -1.0 + phiOffset_)) {
      float phi = rh1.dims[1] + 2 * M_PI;
      KDTreeNodeInfo<size_t, 2> rh3(irechit, float(rechitTable_.get<PF::rechit::Eta>(irechit)), phi);
      eltList.push_back(rh3);
    }
  }

  // Here we define the upper/lower bounds of the 2D space (eta/phi).
  float phimin = -1.0 * M_PI - phiOffset_;
  float phimax = M_PI + phiOffset_;

  // etamin-etamax, phimin-phimax
  KDTreeBox region(-3.0f, 3.0f, phimin, phimax);

  // We may now build the KDTree
  tree_.build(eltList, region);
}

void KDTreeLinkerTrackEcal::searchLinks(const PFTables& pftables) {
  // Most of the code has been taken from LinkByRecHit.cc

  size_t itrack = 0;

  // We iterate over the tracks.
  for (BlockEltSet::iterator it = targetSet_.begin(); it != targetSet_.end(); it++) {
    // We set the multilinks flag of the track to true. It will allow us to
    // use in an optimized way our algo results in the recursive linking algo.
    (*it)->setIsValidMultilinks(true);

    const auto trackPt = trackTable_.get<col::PF::track::Pt>(itrack);
    const auto tracketa = trackTable_.get<col::PF::track::Eta>(itrack);
    const auto trackphi = trackTable_.get<col::PF::track::Phi>(itrack);
    const auto trackx = trackTable_.get<col::PF::track::Posx>(itrack);
    const auto tracky = trackTable_.get<col::PF::track::Posy>(itrack);
    const auto trackz = trackTable_.get<col::PF::track::Posz>(itrack);

    // Estimate the maximal envelope in phi/eta that will be used to find rechit candidates.
    // Same envelope for cap et barrel rechits.
    float range = cristalPhiEtaMaxSize_ * (2.0 + 1.0 / std::min(1., trackPt / 2.));

    // We search for all candidate recHits, ie all recHits contained in the maximal size envelope.
    std::vector<size_t> recHits;
    KDTreeBox trackBox(tracketa - range, tracketa + range, trackphi - range, trackphi + range);
    tree_.search(trackBox, recHits);

    // Here we check all rechit candidates using the non-approximated method.
    for (size_t irecHit : recHits) {
      double rhsizeeta = std::abs(rechitTable_.get<PF::rechit::Corner3eta>(irecHit) -
                                  rechitTable_.get<PF::rechit::Corner1eta>(irecHit));
      double rhsizephi = std::abs(rechitTable_.get<PF::rechit::Corner3phi>(irecHit) -
                                  rechitTable_.get<PF::rechit::Corner1phi>(irecHit));
      if (rhsizephi > M_PI)
        rhsizephi = 2. * M_PI - rhsizephi;

      double deta = std::abs(rechitTable_.get<PF::rechit::Eta>(irecHit) - tracketa);
      double dphi = std::abs(rechitTable_.get<PF::rechit::Phi>(irecHit) - trackphi);
      if (dphi > M_PI)
        dphi = 2. * M_PI - dphi;

      // Find all clusters associated to given rechit
      const auto& rechit_clusters = rechit2ClusterLinksIdx_[irecHit];

      for (const size_t clusteridx : rechit_clusters) {
        //size_t clusteridx = std::distance(fieldClusterVec_.begin(), std::find(fieldClusterVec_.begin(), fieldClusterVec_.end(), *clusterIt));
        auto* clusterPtr = fieldClusterVec_[clusteridx];
        double clusterz = clusterTable_.get<PF::cluster::Posz>(clusteridx);
        int fracsNbr = clusterTable_.get<PF::cluster::fracsNbr>(clusteridx);

        if (clusterTable_.get<PF::cluster::layer>(clusteridx) == PFLayer::ECAL_BARREL) {  // BARREL
          // Check if the track is in the barrel
          if (std::abs(trackz) > 300.)
            continue;

          double _rhsizeeta = rhsizeeta * (2.00 + 1.0 / (fracsNbr * std::min(1., trackPt / 2.)));
          double _rhsizephi = rhsizephi * (2.00 + 1.0 / (fracsNbr * std::min(1., trackPt / 2.)));

          // Check if the track and the cluster are linked
          if (deta < (_rhsizeeta / 2.) && dphi < (_rhsizephi / 2.))
            target2ClusterLinks_[*it].insert(clusterPtr);

        } else {  // ENDCAP

          // Check if the track is in the cap
          if (std::abs(trackz) < 300.)
            continue;
          if (trackz * clusterz < 0.)
            continue;

          double x[5];
          double y[5];
          const double rechit_corner_posx[4] = {rechitTable_.get<PF::rechit::Corner0x>(irecHit),
                                                rechitTable_.get<PF::rechit::Corner1x>(irecHit),
                                                rechitTable_.get<PF::rechit::Corner2x>(irecHit),
                                                rechitTable_.get<PF::rechit::Corner3x>(irecHit)};
          const double rechit_corner_posy[4] = {rechitTable_.get<PF::rechit::Corner0y>(irecHit),
                                                rechitTable_.get<PF::rechit::Corner1y>(irecHit),
                                                rechitTable_.get<PF::rechit::Corner2y>(irecHit),
                                                rechitTable_.get<PF::rechit::Corner3y>(irecHit)};

          for (unsigned jc = 0; jc < 4; ++jc) {
            x[3 - jc] =
                rechit_corner_posx[jc] + (rechit_corner_posx[jc] - rechitTable_.get<PF::rechit::Posx>(irecHit)) *
                                             (1.00 + 0.50 / fracsNbr / std::min(1., trackPt / 2.));
            y[3 - jc] =
                rechit_corner_posy[jc] + (rechit_corner_posy[jc] - rechitTable_.get<PF::rechit::Posy>(irecHit)) *
                                             (1.00 + 0.50 / fracsNbr / std::min(1., trackPt / 2.));
          }

          x[4] = x[0];
          y[4] = y[0];

          bool isinside = TMath::IsInside(trackx, tracky, 5, x, y);

          // Check if the track and the cluster are linked
          if (isinside)
            target2ClusterLinks_[*it].insert(clusterPtr);
        }
      }
    }
    itrack++;
  }
}

void KDTreeLinkerTrackEcal::updatePFBlockEltWithLinks(const PFTables& pftables) {
  //TODO YG : Check if cluster positionREP() is valid ?

  // Here we save in each track the list of phi/eta values of linked clusters.
  for (const auto& track_clusters : target2ClusterLinks_) {
    reco::PFBlockElement* track = track_clusters.first;
    const BlockEltSet& clusters = track_clusters.second;
    reco::PFMultiLinksTC multitracks(true);

    for (reco::PFBlockElement* cluster : clusters) {
      double clusterphi = cluster->clusterRef()->positionREP().phi();
      double clustereta = cluster->clusterRef()->positionREP().eta();

      multitracks.linkedClusters.push_back(std::make_pair(clusterphi, clustereta));
    }

    track->setMultilinks(multitracks);
  }
}

void KDTreeLinkerTrackEcal::clear() {
  targetSet_.clear();
  fieldClusterSet_.clear();
  fieldClusterVec_.clear();

  rechitsSet_.clear();
  rechitsVec_.clear();

  rechit2ClusterLinks_.clear();
  rechit2ClusterLinksIdx_.clear();
  target2ClusterLinks_.clear();

  tree_.clear();
}
