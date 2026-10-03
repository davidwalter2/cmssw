#ifndef Analysis_HitAnalyzer_ParticleProperties_h
#define Analysis_HitAnalyzer_ParticleProperties_h

#include <string>

namespace ana_hitanalyzer {

// PDG mass + mass uncertainty for a charged particle. Both in GeV.
// Used by the V0 / D* candidate producers and the CVH ntuplizers so they
// share a single source of truth for daughter-particle properties.
struct ParticleProperties {
  double mass;
  double massErr;
};

// Look up PDG mass + uncertainty by particle base name. Recognised bases:
//   "mu", "pi", "kaon", "proton", "e".
// Returns {0., 0.} for unknown / empty names; callers are expected to
// guard against this (typically by requiring a non-empty cfi parameter).
ParticleProperties getParticleProperties(const std::string &base);

// Build the Geant4 particle-name string from a base name + track charge.
//   Empty `base` returns empty string -> caller falls back to whatever
//                                        default the propagator was
//                                        constructed with ("mu").
//   "proton"   -> "proton"  (charge > 0) / "anti_proton" (charge <= 0).
//   any other  -> "<base>+" / "<base>-".
std::string g4ParticleName(const std::string &base, int charge);

}  // namespace ana_hitanalyzer

#endif
