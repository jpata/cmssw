#ifndef RecoParticleFlow_PFProducer_PFBlockAlgo_h
#define RecoParticleFlow_PFProducer_PFBlockAlgo_h

#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "RecoParticleFlow/PFProducer/interface/BlockElementImporterBase.h"
#include "RecoParticleFlow/PFProducer/interface/BlockElementLinkerBase.h"
#include "RecoParticleFlow/PFProducer/interface/KDTreeLinkerBase.h"
#include "DataFormats/Common/interface/OwnVector.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace std {
  template <>
  struct hash<std::pair<unsigned int, unsigned int>> {
    typedef std::pair<unsigned int, unsigned int> arg_type;
    typedef unsigned int value_type;
    value_type operator()(const arg_type& arg) const { return arg.first ^ (arg.second << 1); }
  };
}  // namespace std

/// \brief Particle Flow Algorithm
/*!
  \author Colin Bernet (rewrite/refactor by L. Gray)
  \date January 2006 (April 2014) 
*/

//Represents a PF Element on a detector layer
//An element has a local density value
class ElementOnLayer {
public:
  //index to the full element array
  int element_index;
  
  //coordinates of the element on this layer
  float eta, phi;

  //which layer is it (ECAL, HCAL, HF, Tracker, ...)
  int layer_id;

  //local density
  float density;

  //distance to the nearest element with a higher density
  float distance_to_higher;

  //output cluster ID within the layer
  int cluster_id;
};

//A tile contains multiple detector elements
class Tile {
public:
  std::vector<int> elements;
};

constexpr size_t TILES_X = 5;
constexpr size_t TILES_Y = 5;
class TileGrid {
public:
  Tile tiles[TILES_X][TILES_Y];

  void fill(std::vector<ElementOnLayer> elements) {
    int iel = 0;
    for (const auto &el : elements) {
      tiles[getEtaBin(el.eta)][getPhiBin(el.phi)].elements.push_back(iel);
      iel += 1;
    }
  };

  size_t getEtaBin(float eta) const { return 0; };
  size_t getPhiBin(float phi) const { return 0; };
  std::array<size_t, 4> searchBoxEtaPhi(float etaMin, float etaMax, float phiMin, float phiMax) const {
    return {{getEtaBin(etaMin), getEtaBin(etaMax), getPhiBin(phiMin), getPhiBin(phiMax)}};  
  };
};


class PFBlockAlgo {
public:
  // the element list should **always** be a list of (smart) pointers
  typedef std::vector<std::unique_ptr<reco::PFBlockElement>> ElementList;
  //for skipping ranges
  typedef std::array<std::pair<unsigned int, unsigned int>, reco::PFBlockElement::kNBETypes> ElementRanges;

  PFBlockAlgo();

  ~PFBlockAlgo();

  void setLinkers(const std::vector<edm::ParameterSet>&);

  void setImporters(const std::vector<edm::ParameterSet>&, edm::ConsumesCollector&);

  // update event setup info of all linkers
  void updateEventSetup(const edm::EventSetup&);

  // run all of the importers and build KDtrees
  void buildElements(const edm::Event&);

  /// build blocks
  reco::PFBlockCollection findBlocks();

  /// sets debug printout flag
  void setDebug(bool debug) { debug_ = debug; }

  //Retrieve the 1D index of the link tester given the types of both elements
  unsigned getIndex(const reco::PFBlockElement* el1, const reco::PFBlockElement* el2) const;

  //CLUE-specific code
  reco::PFBlockCollection findBlocksCLUE() const;

  //Separate the individual elements per detector layer
  std::vector<std::vector<ElementOnLayer>> buildLayers(const ElementList& elements_) const;

  //Builds the tiles from the elements on one layer
  TileGrid buildTileGrid(const std::vector<ElementOnLayer>& layer) const;

  //Updates the local density of all the elements in the layer
  void calculateLocalDensity(std::vector<ElementOnLayer>& elements, const TileGrid& tiles) const;

  //Updates the distance value to the nearest element with a higher density
  void calculateDistanceToHigher(std::vector<ElementOnLayer>& elements, const TileGrid& tiles) const;

private:
  /// compute missing links in the blocks
  /// (the recursive procedure does not build all links)
  void packLinks(reco::PFBlock& block,
                 const std::unordered_map<std::pair<unsigned int, unsigned int>, double>& links) const;

  /// check whether 2 elements are linked. Returns distance
  inline void link(const reco::PFBlockElement* el1, const reco::PFBlockElement* el2, double& dist) const;

  // the test elements will be transferred to the blocks
  ElementList elements_;
  ElementRanges ranges_;

  /// if true, debug printouts activated
  bool debug_;

  friend std::ostream& operator<<(std::ostream&, const PFBlockAlgo&);
  bool useHO_;

  std::vector<std::unique_ptr<BlockElementImporterBase>> importers_;

  const std::unordered_map<std::string, reco::PFBlockElement::Type> elementTypes_;
  std::vector<std::unique_ptr<BlockElementLinkerBase>> linkTests_;
  unsigned int linkTestSquare_[reco::PFBlockElement::kNBETypes][reco::PFBlockElement::kNBETypes];

  std::vector<std::unique_ptr<KDTreeLinkerBase>> kdtrees_;
};

#endif
