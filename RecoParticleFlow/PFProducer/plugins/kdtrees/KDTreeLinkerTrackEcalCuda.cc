#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerBase.h"
#include "CommonTools/RecoAlgos/interface/KDTreeLinkerAlgo.h"
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"

#include "TMath.h"

#include <cuda_runtime.h>

void track_cluster_dr(
  const float* track_pt,
  const float* track_eta,
  const float* track_phi,
  size_t num_tracks,
  const float* rechit_eta,
  const float* rechit_phi,
  const int* rechit_clusteridx,
  size_t num_rechits,
  int* out_track_clusteridx);

// This class is used to find all links between Tracks and ECAL clusters
// using a KDTree algorithm.
// It is used in PFBlockAlgo.cc in the function links().
class KDTreeLinkerTrackEcalCuda : public KDTreeLinkerBase {
public:
  KDTreeLinkerTrackEcalCuda(const edm::ParameterSet &conf);
  ~KDTreeLinkerTrackEcalCuda() override;

  // With this method, we create the list of psCluster that we want to link.
  void insertTargetElt(reco::PFBlockElement *track) override;

  // Here, we create the list of ecalCluster that we want to link. From ecalCluster
  // and fraction, we will create a second list of rechits that will be used to
  // build the KDTree.
  void insertFieldClusterElt(reco::PFBlockElement *ecalCluster) override;

  // The KDTree building from rechits list.
  void buildTree() override;

  // Here we will iterate over all tracks. For each track intersection point with ECAL,
  // we will search the closest rechits in the KDTree, from rechits we will find the
  // ecalClusters and after that we will check the links between the track and
  // all closest ecalClusters.
  void searchLinks() override;

  // Here, we will store all PS/ECAL founded links in the PFBlockElement class
  // of each psCluster in the PFmultilinks field.
  void updatePFBlockEltWithLinks() override;

  // Here we free all allocated structures.
  void clear() override;

private:

  std::vector<reco::PFBlockElement*> tracks_;
  std::vector<reco::PFBlockElement*> clusters_;

  BlockElt2BlockEltMap tracks_to_clusters_;

  //GPU memory pointers
  float* tracks_pt_;
  float* tracks_eta_;
  float* tracks_phi_;
  float* tracks_x_;
  float* tracks_y_;
  float* tracks_z_;
  size_t num_tracks_;
 
  float* rechits_eta_;
  float* rechits_phi_;
  int* rechits_clusteridx_;
  size_t num_rechits_;
};

// the text name is different so that we can easily
// construct it when calling the factory
DEFINE_EDM_PLUGIN(KDTreeLinkerFactory, KDTreeLinkerTrackEcalCuda, "KDTreeTrackAndECALCudaLinker");

KDTreeLinkerTrackEcalCuda::KDTreeLinkerTrackEcalCuda(const edm::ParameterSet &conf) : KDTreeLinkerBase(conf) {
  std::cout << "KDTreeLinkerTrackEcalCuda constructor" << std::endl;
  cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 100000000);
}

KDTreeLinkerTrackEcalCuda::~KDTreeLinkerTrackEcalCuda() { clear(); }

void KDTreeLinkerTrackEcalCuda::insertTargetElt(reco::PFBlockElement *track) {
  if (track->trackRefPF()->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax).isValid()) {
    tracks_.push_back(track);
  }
}

void KDTreeLinkerTrackEcalCuda::insertFieldClusterElt(reco::PFBlockElement *ecalCluster) {
  clusters_.push_back(ecalCluster);
}

void KDTreeLinkerTrackEcalCuda::buildTree() {
  std::cout << "KDTreeLinkerTrackEcalCuda::buildTree" << std::endl;
  num_tracks_ = tracks_.size();
  std::cout << "num_tracks=" << num_tracks_ << std::endl;
 
  std::vector<float> g_tracks_pt_;
  std::vector<float> g_tracks_eta_;
  std::vector<float> g_tracks_phi_;
  std::vector<float> g_tracks_x_;
  std::vector<float> g_tracks_y_;
  std::vector<float> g_tracks_z_;
  g_tracks_pt_.reserve(num_tracks_);
  g_tracks_eta_.reserve(num_tracks_);
  g_tracks_phi_.reserve(num_tracks_);
  g_tracks_x_.reserve(num_tracks_);
  g_tracks_y_.reserve(num_tracks_);
  g_tracks_z_.reserve(num_tracks_);

  for (const auto* track : tracks_) {
    const auto& trackref = track->trackRefPF();
    const auto& atECAL = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
    const auto& atVertex = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach);
    const auto& posECAL = atECAL.positionREP();

    const float pt = std::sqrt(atVertex.momentum().Vect().Perp2());
    const float eta = posECAL.eta();
    const float phi = posECAL.phi();

    const float x = posECAL.X();
    const float y = posECAL.Y();
    const float z = posECAL.Z();

    g_tracks_pt_.push_back(pt);
    g_tracks_eta_.push_back(eta);
    g_tracks_phi_.push_back(phi);
    g_tracks_x_.push_back(x);
    g_tracks_y_.push_back(y);
    g_tracks_z_.push_back(z);
  }

  cudaCheck(cudaMalloc(&tracks_pt_, sizeof(float)*num_tracks_));
  cudaCheck(cudaMalloc(&tracks_eta_, sizeof(float)*num_tracks_));
  cudaCheck(cudaMalloc(&tracks_phi_, sizeof(float)*num_tracks_));
  cudaCheck(cudaMalloc(&tracks_x_, sizeof(float)*num_tracks_));
  cudaCheck(cudaMalloc(&tracks_y_, sizeof(float)*num_tracks_));
  cudaCheck(cudaMalloc(&tracks_z_, sizeof(float)*num_tracks_));
  cudaCheck(cudaMemcpy(tracks_pt_, g_tracks_pt_.data(), sizeof(float)*g_tracks_pt_.size(), cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(tracks_eta_, g_tracks_eta_.data(), sizeof(float)*g_tracks_eta_.size(), cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(tracks_phi_, g_tracks_phi_.data(), sizeof(float)*g_tracks_phi_.size(), cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(tracks_x_, g_tracks_x_.data(), sizeof(float)*g_tracks_x_.size(), cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(tracks_y_, g_tracks_y_.data(), sizeof(float)*g_tracks_y_.size(), cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(tracks_z_, g_tracks_z_.data(), sizeof(float)*g_tracks_z_.size(), cudaMemcpyDefault));
  cudaCheck(cudaDeviceSynchronize());
  

  std::vector<float> g_rechits_eta;
  std::vector<float> g_rechits_phi;
  std::vector<int> g_rechits_clusteridx;
  
  int icluster = 0;
  num_rechits_ = 0;
  for (const auto* cluster : clusters_) {
    int irechit = 0;
    const auto& fractions = cluster->clusterRef()->recHitFractions();
    for (const auto& frac : fractions) {
      const reco::PFRecHitRef &rh = frac.recHitRef();
      double fract = frac.fraction();

      if ((rh.isNull()) || (fract < cutOffFrac))
        continue;

      const auto& posrep = rh->positionREP();

      g_rechits_eta.push_back(posrep.eta());
      g_rechits_phi.push_back(posrep.phi());
      g_rechits_clusteridx.push_back(icluster);
      std::cout << "rh=" << irechit << " " << posrep.eta() << " " << posrep.phi() << " " << icluster << std::endl;
      num_rechits_++;
      irechit++;
    }
    icluster++;
  }

  std::cout << "num_rechits=" << num_rechits_ << std::endl;
  cudaCheck(cudaMalloc(&rechits_eta_, sizeof(float)*num_rechits_));
  cudaCheck(cudaMalloc(&rechits_phi_, sizeof(float)*num_rechits_));
  cudaCheck(cudaMalloc(&rechits_clusteridx_, sizeof(int)*num_rechits_));
  cudaCheck(cudaMemcpy(rechits_eta_, g_rechits_eta.data(), sizeof(float)*num_rechits_, cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(rechits_phi_, g_rechits_phi.data(), sizeof(float)*num_rechits_, cudaMemcpyDefault));
  cudaCheck(cudaMemcpy(rechits_clusteridx_, g_rechits_clusteridx.data(), sizeof(int)*num_rechits_, cudaMemcpyDefault));
}

void KDTreeLinkerTrackEcalCuda::searchLinks() {

  int* out_track_clusteridx;
  cudaCheck(cudaMalloc(&out_track_clusteridx, sizeof(int)*num_tracks_));
  track_cluster_dr(tracks_pt_, tracks_eta_, tracks_phi_, num_tracks_, rechits_eta_, rechits_phi_, rechits_clusteridx_, num_rechits_, out_track_clusteridx);
  auto track_clusteridx = std::make_unique<int[]>(num_tracks_);
  cudaCheck(cudaMemcpy(track_clusteridx.get(), out_track_clusteridx, sizeof(int)*num_tracks_, cudaMemcpyDeviceToHost));
  cudaCheck(cudaFree(out_track_clusteridx));

  for (int i=0; i<static_cast<int>(num_tracks_); i++) {
    std::cout << track_clusteridx[i] << " "; 
  }
  std::cout << std::endl;

  //this is created by the cuda kernel
  std::vector<std::pair<int, int>> inds;
  inds.push_back(std::make_pair<int, int>(0,0));
  for (const auto ind : inds) {
    const auto ind_track = ind.first;
    const auto ind_cluster = ind.second;
    tracks_to_clusters_[tracks_[ind_track]].insert(clusters_[ind_cluster]);
  }
}

void KDTreeLinkerTrackEcalCuda::updatePFBlockEltWithLinks() {
  for (const auto& track_clusters : tracks_to_clusters_) {
    auto* track = track_clusters.first;
    const auto& clusters = track_clusters.second;

    reco::PFMultiLinksTC multitracks(true);

    for (const auto& cluster : clusters) {
      double clusterphi = cluster->clusterRef()->positionREP().phi();
      double clustereta = cluster->clusterRef()->positionREP().eta();

      multitracks.linkedClusters.push_back(std::make_pair(clusterphi, clustereta));
    }

    track->setMultilinks(multitracks);
  }
}

void KDTreeLinkerTrackEcalCuda::clear() {
  if (tracks_.size() > 0) {
    cudaCheck(cudaFree(tracks_pt_));
    cudaCheck(cudaFree(tracks_eta_));
    cudaCheck(cudaFree(tracks_phi_));
    cudaCheck(cudaFree(tracks_x_));
    cudaCheck(cudaFree(tracks_y_));
    cudaCheck(cudaFree(tracks_z_));
    tracks_.clear();
    cudaCheck(cudaFree(rechits_eta_));
    cudaCheck(cudaFree(rechits_phi_));
    cudaCheck(cudaFree(rechits_clusteridx_));
  }
}
