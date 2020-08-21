#include "RecoParticleFlow/PFProducer/interface/Tables.h"

using namespace edm::soa::col;
namespace edm::soa::col {
TrackTable makeTrackTable(std::set<reco::PFBlockElement*>& targetSet) {
  std::vector<double> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  for (const reco::PFBlockElement* pfelement_track : targetSet) {
    const reco::PFRecTrackRef& trackref = pfelement_track->trackRefPF();

    const reco::PFTrajectoryPoint& atECAL = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ECALShowerMax);
    const reco::PFTrajectoryPoint& atVertex = trackref->extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach);

    pt.push_back(sqrt(atVertex.momentum().Vect().Perp2()));
    eta.push_back(atECAL.positionREP().eta());
    phi.push_back(atECAL.positionREP().phi());
    x.push_back(atECAL.position().X());
    y.push_back(atECAL.position().Y());
    z.push_back(atECAL.position().Z());
  }

  return TrackTable(pt, eta, phi, x, y, z);
}

RecHitTable makeRecHitTable(std::vector<const reco::PFRecHit*> const& objects) {
  return {objects,
          edm::soa::column_fillers(
              PF::rechit::Eta::filler([](reco::PFRecHit const* x) { return x->positionREP().eta(); }),
              PF::rechit::Phi::filler([](reco::PFRecHit const* x) { return x->positionREP().phi(); }),
              PF::rechit::Posx::filler([](reco::PFRecHit const* x) { return x->position().x(); }),
              PF::rechit::Posy::filler([](reco::PFRecHit const* x) { return x->position().y(); }),
              PF::rechit::Posz::filler([](reco::PFRecHit const* x) { return x->position().z(); }),

              PF::rechit::Corner0x::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[0].x(); }),
              PF::rechit::Corner0y::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[0].y(); }),
              PF::rechit::Corner0z::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[0].z(); }),
              PF::rechit::Corner1x::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[1].x(); }),
              PF::rechit::Corner1y::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[1].y(); }),
              PF::rechit::Corner1z::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[1].z(); }),
              PF::rechit::Corner2x::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[2].x(); }),
              PF::rechit::Corner2y::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[2].y(); }),
              PF::rechit::Corner2z::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[2].z(); }),
              PF::rechit::Corner3x::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[3].x(); }),
              PF::rechit::Corner3y::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[3].y(); }),
              PF::rechit::Corner3z::filler([](reco::PFRecHit const* x) { return x->getCornersXYZ()[3].z(); }),

              PF::rechit::Corner0eta::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[0].eta(); }),
              PF::rechit::Corner0phi::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[0].phi(); }),
              PF::rechit::Corner1eta::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[1].eta(); }),
              PF::rechit::Corner1phi::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[1].phi(); }),
              PF::rechit::Corner2eta::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[2].eta(); }),
              PF::rechit::Corner2phi::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[2].phi(); }),
              PF::rechit::Corner3eta::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[3].eta(); }),
              PF::rechit::Corner3phi::filler([](reco::PFRecHit const* x) { return x->getCornersREP()[3].phi(); }))};
}

ClusterTable makeClusterTable(std::vector<reco::PFBlockElement*> const& objects) {
  return {objects,
          edm::soa::column_fillers(
              PF::cluster::Posz::filler([](reco::PFBlockElement* x) { return x->clusterRef()->position().z(); }),
              PF::cluster::fracsNbr::filler(
                  [](reco::PFBlockElement* x) { return x->clusterRef()->recHitFractions().size(); }),
              PF::cluster::layer::filler([](reco::PFBlockElement* x) { return x->clusterRef()->layer(); }))};
}
}