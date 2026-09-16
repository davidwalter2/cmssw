// FlattenedTrackValueMapVectorTableProducer
//
// A reco::Track-keyed instantiation of the stock NanoAOD
// FlattenedValueMapVectorTableProducer template.
//
// PhysicsTools/NanoAOD already provides the reco::Candidate instantiation
// (FlattenedCandValueMapVectorTableProducer), which muons_cff.py uses to store
// exactly this payload on the Muon table -- cvhmergedGlobalIdxs / cvhJacRef /
// cvhMomCov. That covers the joint N-body maker, whose ValueMaps are keyed to
// the VertexCompositeCandidate collection.
//
// It does NOT cover the single-track makers. Their jacobian maps are keyed to
// the input reco::TrackCollection (the B+ bachelor and J/psi-leg makers run on
// bare track collections with no muon association), and reco::Track is not a
// reco::Candidate, so the stock typedef cannot consume them. This adds the one
// missing instantiation rather than duplicating the template.

#include "PhysicsTools/NanoAOD/interface/FlattenedValueMapVectorTableProducer.h"

#include "DataFormats/TrackReco/interface/Track.h"

typedef FlattenedValueMapVectorTableProducer<reco::Track>
    FlattenedTrackValueMapVectorTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FlattenedTrackValueMapVectorTableProducer);
