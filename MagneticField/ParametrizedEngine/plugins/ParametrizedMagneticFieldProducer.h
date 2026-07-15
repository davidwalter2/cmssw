#ifndef ParametrizedMagneticFieldProducer_h
#define ParametrizedMagneticFieldProducer_h

/** \class ParametrizedMagneticFieldProducer
 *
 *   Description: Producer for the Parametrized Magnetic Field
 *
 *  \author Massimiliano Chiorboli, updated NA 03/08
 */

#include "FWCore/Framework/interface/ESProducer.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <memory>

class IdealMagneticFieldRecord;

namespace magneticfield {
  class ParametrizedMagneticFieldProducer : public edm::ESProducer {
  public:
    ParametrizedMagneticFieldProducer(const edm::ParameterSet&);
    ~ParametrizedMagneticFieldProducer() override;

    std::shared_ptr<MagneticField> produce(const IdealMagneticFieldRecord&);
    edm::ParameterSet pset;

  private:
    // The parametrized field depends only on the configuration PSet (version
    // + parameter block), never on run conditions -- unlike Auto*Producer it
    // does not follow the magnet current. Build it once and hand the SAME
    // object to every IOV of IdealMagneticFieldRecord (which rolls every run
    // through its RunInfoRcd dependency). This keeps the product pointer
    // stable across run boundaries, which the CVH G4 world relies on: it
    // captures the raw MagneticField* at job start and cannot be re-pointed.
    std::shared_ptr<MagneticField> field_;
  };
}  // namespace magneticfield

#endif
