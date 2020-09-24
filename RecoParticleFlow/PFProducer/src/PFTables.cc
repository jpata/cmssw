#include "RecoParticleFlow/PFProducer/interface/PFTables.h"

using namespace edm::soa::col;
namespace edm::soa {

  TrackTableVertex makeTrackTableVertex(std::vector<reco::PFBlockElement*> const& objects) {
    return {objects,
            edm::soa::column_fillers(
                pf::track::ExtrapolationValid::filler([](reco::PFBlockElement* x) {
                  return x->trackRefPF()->extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach).isValid();
                }),
                pf::track::Pt::filler([](reco::PFBlockElement* x) {
                  return sqrt(x->trackRefPF()
                                  ->extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach)
                                  .momentum()
                                  .Vect()
                                  .Perp2());
                }),
                pf::track::IsLinkedToDisplacedVertex::filler(
                    [](reco::PFBlockElement* x) { return x->isLinkedToDisplacedVertex(); }),
                pf::track::TrackType_FROM_GAMMACONV::filler(
                    [](reco::PFBlockElement* x) { return x->trackType(reco::PFBlockElement::T_FROM_GAMMACONV); }),
                pf::track::TrackType_FROM_V0::filler(
                    [](reco::PFBlockElement* x) { return x->trackType(reco::PFBlockElement::T_FROM_V0); }),
                pf::track::V0RefIsNonNull::filler([](reco::PFBlockElement* x) { return x->V0Ref().isNonnull(); }),
                pf::track::V0RefKey::filler([](reco::PFBlockElement* x) { return refToElementID(x->V0Ref()); }),
                pf::track::KfTrackRefIsNonNull::filler(
                    [](reco::PFBlockElement* x) { return x->trackRefPF()->trackRef().isNonnull(); }),
                pf::track::KfTrackRefKey::filler(
                    [](reco::PFBlockElement* x) { return refToElementID(x->trackRefPF()->trackRef()); }),
                pf::track::KfTrackRefBaseKey::filler([](reco::PFBlockElement* x) {
                  return refToElementID(reco::TrackBaseRef(x->trackRefPF()->trackRef()));
                }),
                pf::track::DisplacedVertexRef_TO_DISP_IsNonNull::filler([](reco::PFBlockElement* x) {
                  return x->displacedVertexRef(reco::PFBlockElement::T_TO_DISP).isNonnull();
                }),
                pf::track::DisplacedVertexRef_TO_DISP_Key::filler([](reco::PFBlockElement* x) {
                  return refToElementID(x->displacedVertexRef(reco::PFBlockElement::T_TO_DISP));
                }),
                pf::track::DisplacedVertexRef_FROM_DISP_IsNonNull::filler([](reco::PFBlockElement* x) {
                  return x->displacedVertexRef(reco::PFBlockElement::T_FROM_DISP).isNonnull();
                }),
                pf::track::DisplacedVertexRef_FROM_DISP_Key::filler([](reco::PFBlockElement* x) {
                  return refToElementID(x->displacedVertexRef(reco::PFBlockElement::T_FROM_DISP));
                }))};
  }

  ConvRefTable makeConvRefTable(const std::vector<reco::ConversionRef>& convrefs) {
    return {convrefs,
            edm::soa::column_fillers(
                pf::track::ConvRefIsNonNull::filler([](const reco::ConversionRef& x) { return x.isNonnull(); }),
                pf::track::ConvRefKey::filler([](const reco::ConversionRef& x) { return refToElementID(x); }))};
  }

  ConvBremTable makeConvBremTable(const std::vector<reco::PFRecTrackRef>& convbrems) {
    return {convbrems,
            edm::soa::column_fillers(
                pf::track::ConvBremRefKey::filler([](reco::PFRecTrackRef x) { return refToElementID(x->trackRef()); }),
                pf::track::ConvBremRefBaseKey::filler(
                    [](reco::PFRecTrackRef x) { return refToElementID(reco::TrackBaseRef(x->trackRef())); }))};
  }

  BremTable makeBremTable(const std::vector<const reco::PFBlockElementBrem*>& brems) {
    return {
        brems,
        edm::soa::column_fillers(
            pf::track::Pt::filler([](const reco::PFBlockElementBrem* x) {
              return sqrt(
                  x->trackPF().extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach).momentum().Vect().Perp2());
            }),
            pf::track::GsfTrackRefPFIsNonNull::filler(
                [](const reco::PFBlockElementBrem* x) { return x->GsftrackRefPF().isNonnull(); }),
            pf::track::GsfTrackRefPFKey::filler(
                [](const reco::PFBlockElementBrem* x) { return refToElementID(x->GsftrackRefPF()); }))};
  }

  TrackTableExtrapolation makeTrackTable(std::vector<reco::PFBlockElement*> const& objects,
                                         reco::PFTrajectoryPoint::LayerType layerType) {
    std::vector<reco::PFRecTrackRef> rectracks;
    rectracks.reserve(objects.size());
    for (const auto& obj : objects) {
      rectracks.push_back(obj->trackRefPF());
    }
    return makeTrackTable(rectracks, layerType);
  }

  //this should be unified with the function above
  TrackTableExtrapolation makeTrackTable(std::vector<const reco::PFBlockElementGsfTrack*> const& objects,
                                         reco::PFTrajectoryPoint::LayerType layerType) {
    std::vector<const reco::GsfPFRecTrack*> rectracks;
    rectracks.reserve(objects.size());
    for (const auto* obj : objects) {
      rectracks.push_back(&(obj->GsftrackPF()));
    }
    return makeTrackTable(rectracks, layerType);
  }

  TrackTableExtrapolation makeTrackTable(std::vector<const reco::PFBlockElementBrem*> const& objects,
                                         reco::PFTrajectoryPoint::LayerType layerType) {
    std::vector<const reco::PFRecTrack*> rectracks;
    rectracks.reserve(objects.size());
    for (const auto* obj : objects) {
      rectracks.push_back(&(obj->trackPF()));
    }
    return makeTrackTable(rectracks, layerType);
  }

  template <class RecTrackType>
  TrackTableExtrapolation makeTrackTable(std::vector<RecTrackType> const& objects,
                                         reco::PFTrajectoryPoint::LayerType layerType) {
    return {objects,
            edm::soa::column_fillers(pf::track::ExtrapolationValid::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).isValid();
                                     }),
                                     pf::track::Eta::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).positionREP().eta();
                                     }),
                                     pf::track::Phi::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).positionREP().phi();
                                     }),
                                     pf::track::Posx::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).positionREP().X();
                                     }),
                                     pf::track::Posy::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).positionREP().Y();
                                     }),
                                     pf::track::Posz::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).positionREP().Z();
                                     }),
                                     pf::track::PosR::filler([layerType](RecTrackType x) {
                                       return x->extrapolatedPoint(layerType).positionREP().R();
                                     }))};
  }

  RecHitTable makeRecHitTable(std::vector<const reco::PFRecHitFraction*> const& objects) {
    return {
        objects,
        edm::soa::column_fillers(
            pf::rechit::DetIdValue::filler([](reco::PFRecHitFraction const* x) { return x->recHitRef()->detId(); }),
            pf::rechit::Fraction::filler([](reco::PFRecHitFraction const* x) { return x->fraction(); }),
            pf::rechit::Eta::filler([](reco::PFRecHitFraction const* x) { return x->recHitRef()->positionREP().eta(); }),
            pf::rechit::Phi::filler([](reco::PFRecHitFraction const* x) { return x->recHitRef()->positionREP().phi(); }),
            pf::rechit::Posx::filler([](reco::PFRecHitFraction const* x) { return x->recHitRef()->position().x(); }),
            pf::rechit::Posy::filler([](reco::PFRecHitFraction const* x) { return x->recHitRef()->position().y(); }),
            pf::rechit::Posz::filler([](reco::PFRecHitFraction const* x) { return x->recHitRef()->position().z(); }),

            pf::rechit::Corner0x::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[0].x(); }),
            pf::rechit::Corner0y::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[0].y(); }),
            pf::rechit::Corner0z::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[0].z(); }),
            pf::rechit::Corner1x::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[1].x(); }),
            pf::rechit::Corner1y::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[1].y(); }),
            pf::rechit::Corner1z::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[1].z(); }),
            pf::rechit::Corner2x::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[2].x(); }),
            pf::rechit::Corner2y::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[2].y(); }),
            pf::rechit::Corner2z::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[2].z(); }),
            pf::rechit::Corner3x::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[3].x(); }),
            pf::rechit::Corner3y::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[3].y(); }),
            pf::rechit::Corner3z::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[3].z(); }),

            pf::rechit::Corner0eta::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[0].eta(); }),
            pf::rechit::Corner0phi::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[0].phi(); }),
            pf::rechit::Corner1eta::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[1].eta(); }),
            pf::rechit::Corner1phi::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[1].phi(); }),
            pf::rechit::Corner2eta::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[2].eta(); }),
            pf::rechit::Corner2phi::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[2].phi(); }),
            pf::rechit::Corner3eta::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[3].eta(); }),
            pf::rechit::Corner3phi::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersREP()[3].phi(); }),

            pf::rechit::Corner0xBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[0].basicVector().x(); }),
            pf::rechit::Corner0yBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[0].basicVector().y(); }),
            pf::rechit::Corner1xBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[1].basicVector().x(); }),
            pf::rechit::Corner1yBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[1].basicVector().y(); }),
            pf::rechit::Corner2xBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[2].basicVector().x(); }),
            pf::rechit::Corner2yBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[2].basicVector().y(); }),
            pf::rechit::Corner3xBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[3].basicVector().x(); }),
            pf::rechit::Corner3yBV::filler(
                [](reco::PFRecHitFraction const* x) { return x->recHitRef()->getCornersXYZ()[3].basicVector().y(); })

                )};
  }

  SuperClusterRecHitTable makeSuperClusterRecHitTable(std::vector<const std::pair<DetId, float>*> const& objects) {
    return {objects,
            edm::soa::column_fillers(
                pf::rechit::DetIdValue::filler([](const std::pair<DetId, float>* x) { return x->first; }),
                pf::rechit::Fraction::filler([](const std::pair<DetId, float>* x) { return x->second; }))};
  }

  ClusterTable makeClusterTable(std::vector<const reco::PFBlockElementCluster*> const& objects) {
    return {objects,
            edm::soa::column_fillers(pf::cluster::Eta::filler([](const reco::PFBlockElementCluster* x) {
                                       return x->clusterRef()->positionREP().eta();
                                     }),
                                     pf::cluster::Phi::filler([](const reco::PFBlockElementCluster* x) {
                                       return x->clusterRef()->positionREP().phi();
                                     }),
                                     pf::cluster::Posx::filler([](const reco::PFBlockElementCluster* x) {
                                       return x->clusterRef()->position().x();
                                     }),
                                     pf::cluster::Posy::filler([](const reco::PFBlockElementCluster* x) {
                                       return x->clusterRef()->position().y();
                                     }),
                                     pf::cluster::Posz::filler([](const reco::PFBlockElementCluster* x) {
                                       return x->clusterRef()->position().z();
                                     }),
                                     pf::cluster::FracsNbr::filler([](const reco::PFBlockElementCluster* x) {
                                       return x->clusterRef()->recHitFractions().size();
                                     }),
                                     pf::cluster::Layer::filler(
                                         [](const reco::PFBlockElementCluster* x) { return x->clusterRef()->layer(); }),
                                     pf::cluster::SCRefKey::filler([](const reco::PFBlockElementCluster* x) {
                                       return refToElementID(x->superClusterRef());
                                     }))};
  }

  SuperClusterTable makeSuperClusterTable(std::vector<const reco::PFBlockElementSuperCluster*> const& objects) {
    return {objects,
            edm::soa::column_fillers(pf::cluster::Eta::filler([](const reco::PFBlockElementSuperCluster* x) {
                                       return x->superClusterRef()->position().eta();
                                     }),
                                     pf::cluster::Phi::filler([](const reco::PFBlockElementSuperCluster* x) {
                                       return x->superClusterRef()->position().phi();
                                     }),
                                     pf::cluster::SCRefKey::filler([](const reco::PFBlockElementSuperCluster* x) {
                                       return refToElementID(x->superClusterRef());
                                     }))};
  }

  GSFTable makeGSFTable(std::vector<const reco::PFBlockElementGsfTrack*> const& objects) {
    return {
        objects,
        edm::soa::column_fillers(
            pf::track::Pt::filler([](const reco::PFBlockElementGsfTrack* x) {
              return sqrt(
                  x->GsftrackPF().extrapolatedPoint(reco::PFTrajectoryPoint::ClosestApproach).momentum().Vect().Perp2());
            }),
            pf::track::GsfTrackRefPFIsNonNull::filler(
                [](const reco::PFBlockElementGsfTrack* x) { return x->GsftrackRefPF().isNonnull(); }),
            pf::track::GsfTrackRefPFKey::filler(
                [](const reco::PFBlockElementGsfTrack* x) { return refToElementID(x->GsftrackRefPF()); }),
            pf::track::KfPFRecTrackRefIsNonNull::filler([](const reco::PFBlockElementGsfTrack* x) {
              return x->GsftrackRefPF()->kfPFRecTrackRef().isNonnull();
            }),
            pf::track::KfPFRecTrackRefKey::filler([](const reco::PFBlockElementGsfTrack* x) {
              return refToElementID(x->GsftrackRefPF()->kfPFRecTrackRef());
            }),

            pf::track::KfTrackRefIsNonNull::filler([](const reco::PFBlockElementGsfTrack* x) {
              const auto& r1 = x->GsftrackRefPF();
              bool ret = false;
              if (r1.isNonnull()) {
                const auto r2 = r1->kfPFRecTrackRef();
                if (r2.isNonnull()) {
                  const auto r3 = r2->trackRef();
                  ret = r3.isNonnull();
                }
              }
              return ret;
            }),
            pf::track::KfTrackRefKey::filler([](const reco::PFBlockElementGsfTrack* x) {
              const auto& r1 = x->GsftrackRefPF();
              edm::ElementID ret;
              if (r1.isNonnull()) {
                const auto r2 = r1->kfPFRecTrackRef();
                if (r2.isNonnull()) {
                  const auto r3 = r2->trackRef();
                  ret = refToElementID(r3);
                }
              }
              return ret;
            }),

            pf::track::TrackType_FROM_GAMMACONV::filler([](const reco::PFBlockElementGsfTrack* x) {
              return x->trackType(reco::PFBlockElement::T_FROM_GAMMACONV);
            }),
            pf::track::GsfTrackRefPFTrackId::filler(
                [](const reco::PFBlockElementGsfTrack* x) { return x->GsftrackRefPF()->trackId(); }))};
  }
}  // namespace edm::soa