// WMass: private copy of the CMSSW_10_6_26 PUPPI code (interface/PuppiCandidate.h) for the PV-robust
// PUPPI + DeepMET of WmassNanoProd_10_6_26. The PV-robust DeepMET models were trained
// on the weights this version produces; the PuppiContainer of the release changed the
// alpha RMS/median inputs and the neutral treatment since, so the release class would
// change the inputs of the model. Only the namespace and the include paths differ.
#ifndef WMASS_PUPPI106_CommonTools_PileupAlgos_PuppiCandidate_h
#define WMASS_PUPPI106_CommonTools_PileupAlgos_PuppiCandidate_h

namespace puppi106 {

struct PuppiCandidate {
  double pt{0};
  double eta{0};
  double phi{0};
  double m{0};
  double rapidity{0};
  double px{0};
  double py{0};
  double pz{0};
  double e{0};
  int id{0};
  int puppi_register{0};
};

}  // namespace puppi106

#endif
