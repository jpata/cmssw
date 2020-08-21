#ifndef RecoParticleFlow_PFProducer_Tables_h
#define RecoParticleFlow_PFProducer_Tables_h
#include <set>
#include <vector>

#include "DataFormats/ParticleFlowReco/interface/PFBlockElement.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "CommonTools/Utils/interface/KinematicTables.h"

namespace edm {
  namespace soa {
    namespace col {
      namespace PF {
        namespace track {
          SOA_DECLARE_COLUMN(Pt, double, "pt");
          SOA_DECLARE_COLUMN(Eta, float, "eta");
          SOA_DECLARE_COLUMN(Phi, float, "phi");
          SOA_DECLARE_COLUMN(Posx, double, "Posx");
          SOA_DECLARE_COLUMN(Posy, double, "Posy");
          SOA_DECLARE_COLUMN(Posz, double, "Posz");
        };  // namespace track
        namespace rechit {
          SOA_DECLARE_COLUMN(Eta, double, "eta");
          SOA_DECLARE_COLUMN(Phi, double, "phi");
          SOA_DECLARE_COLUMN(Posx, double, "Posx");
          SOA_DECLARE_COLUMN(Posy, double, "Posy");
          SOA_DECLARE_COLUMN(Posz, double, "Posz");

          SOA_DECLARE_COLUMN(Corner0x, double, "Corner0x");
          SOA_DECLARE_COLUMN(Corner0y, double, "Corner0y");
          SOA_DECLARE_COLUMN(Corner0z, double, "Corner0z");
          SOA_DECLARE_COLUMN(Corner1x, double, "Corner1x");
          SOA_DECLARE_COLUMN(Corner1y, double, "Corner1y");
          SOA_DECLARE_COLUMN(Corner1z, double, "Corner1z");
          SOA_DECLARE_COLUMN(Corner2x, double, "Corner2x");
          SOA_DECLARE_COLUMN(Corner2y, double, "Corner2y");
          SOA_DECLARE_COLUMN(Corner2z, double, "Corner2z");
          SOA_DECLARE_COLUMN(Corner3x, double, "Corner3x");
          SOA_DECLARE_COLUMN(Corner3y, double, "Corner3y");
          SOA_DECLARE_COLUMN(Corner3z, double, "Corner3z");

          SOA_DECLARE_COLUMN(Corner0eta, double, "Corner0eta");
          SOA_DECLARE_COLUMN(Corner0phi, double, "Corner0phi");
          SOA_DECLARE_COLUMN(Corner1eta, double, "Corner1eta");
          SOA_DECLARE_COLUMN(Corner1phi, double, "Corner1phi");
          SOA_DECLARE_COLUMN(Corner2eta, double, "Corner2eta");
          SOA_DECLARE_COLUMN(Corner2phi, double, "Corner2phi");
          SOA_DECLARE_COLUMN(Corner3eta, double, "Corner3eta");
          SOA_DECLARE_COLUMN(Corner3phi, double, "Corner3phi");
        }  // namespace rechit
        namespace cluster {
          SOA_DECLARE_COLUMN(Posz, double, "Posz");
          SOA_DECLARE_COLUMN(fracsNbr, int, "fracsNbr");
          SOA_DECLARE_COLUMN(layer, PFLayer::Layer, "layer");
        }  // namespace cluster
      }    // namespace PF

      using TrackTable =
          Table<PF::track::Pt, PF::track::Eta, PF::track::Phi, PF::track::Posx, PF::track::Posy, PF::track::Posz>;
      using TrackTableView = TableView<PF::track::Pt, PF::track::Eta, PF::track::Phi, PF::track::Posx, PF::track::Posy, PF::track::Posz>;

      using RecHitTable = Table<PF::rechit::Eta,
                                PF::rechit::Phi,
                                PF::rechit::Posx,
                                PF::rechit::Posy,
                                PF::rechit::Posz,
                                PF::rechit::Corner0x,
                                PF::rechit::Corner0y,
                                PF::rechit::Corner0z,
                                PF::rechit::Corner1x,
                                PF::rechit::Corner1y,
                                PF::rechit::Corner1z,
                                PF::rechit::Corner2x,
                                PF::rechit::Corner2y,
                                PF::rechit::Corner2z,
                                PF::rechit::Corner3x,
                                PF::rechit::Corner3y,
                                PF::rechit::Corner3z,
                                PF::rechit::Corner0eta,
                                PF::rechit::Corner0phi,
                                PF::rechit::Corner1eta,
                                PF::rechit::Corner1phi,
                                PF::rechit::Corner2eta,
                                PF::rechit::Corner2phi,
                                PF::rechit::Corner3eta,
                                PF::rechit::Corner3phi>;
      using ClusterTable = Table<PF::cluster::Posz, PF::cluster::fracsNbr, PF::cluster::layer>;

      extern TrackTable makeTrackTable(std::set<reco::PFBlockElement*>& targetSet);
      extern RecHitTable makeRecHitTable(std::vector<const reco::PFRecHit*> const& objects);
      extern ClusterTable makeClusterTable(std::vector<reco::PFBlockElement*> const& objects);
    }      // namespace col
  }        // namespace soa
}  // namespace edm
#endif