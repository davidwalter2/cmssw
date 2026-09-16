// WMass: private copy of the CMSSW_10_6_26 PUPPI code (interface/RecoObj.h) for the PV-robust
// PUPPI + DeepMET of WmassNanoProd_10_6_26. The PV-robust DeepMET models were trained
// on the weights this version produces; the PuppiContainer of the release changed the
// alpha RMS/median inputs and the neutral treatment since, so the release class would
// change the inputs of the model. Only the namespace and the include paths differ.
#ifndef WMASS_PUPPI106_CommonTools_PileupAlgos_PUPPI_RECOOBJ_HH
#define WMASS_PUPPI106_CommonTools_PileupAlgos_PUPPI_RECOOBJ_HH

namespace puppi106 {

class RecoObj 
{
public:
      RecoObj():
	pt(0), eta(0), phi(0), m(0),
	id(0),pfType(-1),vtxId(-1),
	trkChi2(0),vtxChi2(0),
	time(0),depth(0),
	expProb(0),expChi2PU(0),expChi2(0),
	dZ(0),d0(0),charge(0)
    {}
    ~RecoObj(){}
    
    float         pt, eta, phi, m, rapidity;  // kinematics
    int           id;
    int           pfType;
    int           vtxId;               // Vertex Id from Vertex Collection
    float         trkChi2;             // Track Chi2
    float         vtxChi2;             // Vertex Chi2
    float         time,depth;    // Usefule Info
    float         expProb;
    float         expChi2PU;
    float         expChi2;
    float         dZ;
    float         d0;
    int           charge;
};
}  // namespace puppi106

#endif
