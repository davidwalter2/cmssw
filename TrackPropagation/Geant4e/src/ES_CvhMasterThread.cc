#include "TrackPropagation/Geant4e/interface/CvhMasterThread.h"
#include "FWCore/Utilities/interface/typelookup.h"

// Register CvhMasterThread as an EventSetup data product so it can be
// produced on CvhMasterRecord and consumed via esConsumes by the CVH makers.
TYPELOOKUP_DATA_REG(CvhMasterThread);
