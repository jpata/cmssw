#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"

#include <iomanip>
#include <sstream>

using namespace std;
using namespace reco;

void PFBlock::addElement(PFBlockElement* element) {
  element->setIndex(elements_.size());
  element->lock();
  elements_.push_back(element->clone());
}

void PFBlock::bookLinkData() {}

void PFBlock::setLink(unsigned i1, unsigned i2, double Dist, LinkData& linkData, LinkTest test) const {
  assert(test < LINKTEST_ALL);

  unsigned index = 0;
  bool ok = matrix2vector(i1, i2, index);

  if (ok) {
    //ignore the  -1, -1 pair
    if (Dist > -0.5) {
      Link& l = linkData[index];
      l.distance = Dist;
      l.test |= (1 << test);
    } else  //delete if existing
    {
      LinkData::iterator it = linkData.find(index);
      if (it != linkData.end())
        linkData.erase(it);
    }

  } else {
    assert(0);
  }
}

// void PFBlock::lock(unsigned i, LinkData& linkData ) const {

//   assert( linkData.size() == linkDataSize() );

//   for(unsigned j=0; j<elements_.size(); j++) {

//     if(i==j) continue;

//     unsigned index = 0;
//     bool ok =  matrix2vector(i,j, index);
//     if(ok)
//       linkData[index] = -1;
//     else
//       assert(0);
//   }
// }

void PFBlock::associatedElements(unsigned i,
                                 const LinkData& linkData,
                                 multimap<double, unsigned>& sortedAssociates,
                                 PFBlockElement::Type type,
                                 LinkTest test) const {
  sortedAssociates.clear();

  // i is too large
  if (i > elements_.size())
    return;
  // assert(i>=0); i >= 0, since i is unsigned

  for (unsigned ie = 0; ie < elements_.size(); ie++) {
    // considered element itself
    if (ie == i) {
      continue;
    }
    // not the right type
    if (type != PFBlockElement::NONE && elements_[ie].type() != type) {
      continue;
    }

    // Order the elements by increasing distance !

    unsigned index = 0;
    if (!matrix2vector(i, ie, index))
      continue;

    double c2 = -1;
    LinkData::const_iterator it = linkData.find(index);
    if (it != linkData.end() && (((1 << test) & it->second.test) != 0 || (test == LINKTEST_ALL)))
      c2 = it->second.distance;

    // not associated
    if (c2 < 0) {
      continue;
    }

    sortedAssociates.insert(pair<double, unsigned>(c2, ie));
  }
}

bool PFBlock::matrix2vector(unsigned iindex, unsigned jindex, unsigned& index) const {
  unsigned size = elements_.size();
  if (iindex == jindex || iindex >= size || jindex >= size) {
    return false;
  }

  if (iindex > jindex)
    std::swap(iindex, jindex);

  index = jindex - iindex - 1;

  if (iindex > 0) {
    index += iindex * size;
    unsigned missing = iindex * (iindex + 1) / 2;
    index -= missing;
  }

  return true;
}

double PFBlock::dist(unsigned ie1, unsigned ie2, const LinkData& linkData) const {
  double Dist = -1;

  unsigned index = 0;
  if (!matrix2vector(ie1, ie2, index))
    return Dist;
  LinkData::const_iterator it = linkData.find(index);
  if (it != linkData.end())
    Dist = it->second.distance;

  return Dist;
}

ostream& reco::operator<<(ostream& out, const reco::PFBlock& block) {
  if (!out)
    return out;
  out << endl << "\"elements\": [" << endl; 
  for (unsigned i = 0; i < block.elements().size(); i++) {
    out << "  " << block.elements()[i] << "," << endl;
  }
  out << "]," << endl;

  out << "\"linkData\": {"; 
  if (!block.linkData().empty()) {
    for (unsigned i = 0; i < block.elements().size(); i++) {
      for (unsigned j = 0; j < block.elements().size(); j++) {
        double Dist = block.dist(i, j, block.linkData());
        if (Dist > -0.5) {
          out << "(" << i << "," << j << "):" << Dist << ", ";
        }
      }
    }
  }
  out << "},";
  
  return out;
}

unsigned PFBlock::linkDataSize() const {
  unsigned n = elements_.size();

  // number of possible undirected links between n elements.
  // reflective links impossible.

  return n * (n - 1) / 2;
}
