#ifndef TrackPropagation_Geant4e_CvhMasterRecord_h
#define TrackPropagation_Geant4e_CvhMasterRecord_h

// EventSetup record carrying the shared CVH Geant4 master (CvhMasterThread).
//
// The master owns a dedicated G4 thread + G4MTRunManagerKernel (a Geant4
// process-global singleton), so exactly ONE may exist per job. Making it an ES
// product on this record -- consumed via esConsumes by every CVH residual
// maker -- lets any number/mix of maker instances (single-track nominal/ideal/
// bs, two-track, V0, B->J/psi K) coexist in one job while sharing that single
// master, instead of each producer building its own via edm::GlobalCache
// (which is per-module-label and therefore collides on the second instance).
//
// It depends on IdealGeometryRecord (DDD/DD4hep world) and
// IdealMagneticFieldRecord (the labelled master field): the CvhMasterESProducer
// consumes those to build the master. Both are single-IOV over a normal job, so
// the master is produced exactly once.

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "FWCore/Utilities/interface/mplVector.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

class CvhMasterRecord
    : public edm::eventsetup::DependentRecordImplementation<
          CvhMasterRecord,
          edm::mpl::Vector<IdealGeometryRecord, IdealMagneticFieldRecord>> {};

#endif
