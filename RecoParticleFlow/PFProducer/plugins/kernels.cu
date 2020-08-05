#include <cuda_runtime.h>
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/SimpleVector.h"

static const float cristalPhiEtaMaxSize_ = 0.04f;

__global__ void kernel_track_cluster_dr(
  const float* track_pt,
  const float* track_eta,
  const float* track_phi,
  size_t num_tracks,
  const float* rechit_eta,
  const float* rechit_phi,
  const int* rechit_clusteridx,
  size_t num_rechits,
  int* out_track_clusteridx)
{
  //process each track in a different CUDA thread
  for (size_t i=blockDim.x*blockIdx.x+threadIdx.x; i<num_tracks; i+=blockDim.x*gridDim.x) {

    const float range = cristalPhiEtaMaxSize_ * (2.0f + 1.0f / std::min(1.0f, track_pt[i] / 2.0f));
    printf("A track=%d pt=%.2f range=%.2f\n", (int)i, track_pt[i], range);

    // find all rechits in box around track with track_eta+-range, track_phi+-range
    cms::cuda::SimpleVector<int> rechit_inds;
    for (size_t j=0; j<num_rechits; j++) {
      const float deta = std::abs(track_eta[i] - rechit_eta[j]);
      const float dphi = std::abs(track_phi[i] - rechit_phi[j]);
      const bool match = (deta<range && dphi<range);
      if (match) {
        rechit_inds.push_back_unsafe(j);
      }
    }

    // loop over clusters associated to each found rechit
    for (const auto rh_idx : rechit_inds) {
    }
  } //loop over tracks
  

}


 
void track_cluster_dr(
  const float* track_pt,
  const float* track_eta,
  const float* track_phi,
  size_t num_tracks,
  const float* rechit_eta,
  const float* rechit_phi,
  const int* rechit_clusteridx,
  size_t num_rechits,
  int* out_track_clusteridx)
{
  kernel_track_cluster_dr<<<1, 1>>>(track_pt, track_eta, track_phi, num_tracks, rechit_eta, rechit_phi, rechit_clusteridx, num_rechits, out_track_clusteridx);
  cudaCheck(cudaDeviceSynchronize());
  cudaCheck(cudaGetLastError());
}

