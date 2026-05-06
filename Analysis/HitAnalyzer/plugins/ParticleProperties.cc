#include "Analysis/HitAnalyzer/interface/ParticleProperties.h"

namespace ana_hitanalyzer {

// PDG 2024 values (https://pdg.lbl.gov)
//   pion   : 139.57039 +/- 0.00018 MeV
//   kaon   : 493.677   +/- 0.013   MeV
//   muon   : 105.6583755 +/- 0.0000023 MeV
//   proton : 938.27208943 +/- 0.00000029 MeV
//   electron: 0.51099895069 +/- 0.00000000016 MeV
ParticleProperties getParticleProperties(const std::string &base) {
  if (base == "mu")     return {0.1056583755,    0.0000000023};
  if (base == "pi")     return {0.13957039,      0.00000018};
  if (base == "kaon")   return {0.493677,        0.000013};
  if (base == "proton") return {0.93827208943,   0.00000000029};
  if (base == "e")      return {0.000510998950,  0.0000000000040};
  return {0.0, 0.0};
}

std::string g4ParticleName(const std::string &base, int charge) {
  if (base.empty()) return std::string();
  if (base == "proton") return (charge > 0) ? "proton" : "anti_proton";
  return base + (charge > 0 ? "+" : "-");
}

}  // namespace ana_hitanalyzer
