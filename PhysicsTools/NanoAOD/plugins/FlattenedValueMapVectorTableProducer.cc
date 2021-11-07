#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"
#include "PhysicsTools/NanoAOD/interface/FlattenedValueMapVectorTableProducer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
typedef SimpleFlatTableProducer<reco::Candidate> SimpleCandidateFlatTableProducer;
typedef FlattenedValueMapVectorTableProducer<reco::Candidate> FlattenedCandValueMapVectorTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(FlattenedCandValueMapVectorTableProducer);
