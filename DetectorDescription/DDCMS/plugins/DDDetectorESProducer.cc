// -*- C++ -*-
//
// Package:    DetectorDescription/DDDetectorESProducer
// Class:      DDDetectorESProducer
//
/**\class DDDetectorESProducer

 Description: Produce Detector description

 Implementation:
     Detector is described in XML
*/
//
// Original Author:  Ianna Osborne
//         Created:  Fri, 07 Dec 2018 11:20:52 GMT
//
//

#include <memory>

#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Concurrency/interface/SharedResourceNames.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/Common/interface/FileBlob.h"
#include "CondFormats/DataRecord/interface/MFGeometryFileRcd.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/GeometryFileRcd.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DetectorDescription/DDCMS/interface/DDDetector.h"
#include "DD4hep/Detector.h"

using namespace std;
using namespace cms;
using namespace edm;

class DDDetectorESProducer : public ESProducer, public EventSetupRecordIntervalFinder {
public:
  DDDetectorESProducer(const ParameterSet&);
  ~DDDetectorESProducer() override;

  using ReturnType = shared_ptr<DDDetector>;
  using Detector = dd4hep::Detector;

  ReturnType produceGeom(const IdealGeometryRecord&);
  ReturnType produceMagField(const IdealMagneticFieldRecord&);
  ReturnType produce();
  static void fillDescriptions(ConfigurationDescriptions&);

protected:
  void setIntervalFor(const eventsetup::EventSetupRecordKey&, const IOVSyncValue&, ValidityInterval&) override;

private:
  const bool fromDB_;
  const string appendToDataLabel_;
  const string confGeomXMLFiles_;
  const string rootDDName_;
  const string label_;
  edm::ESGetToken<FileBlob, MFGeometryFileRcd> mfToken_;
  edm::ESGetToken<FileBlob, GeometryFileRcd> geomToken_;
  // The XML-driven detector is a pure function of (appendToDataLabel, confGeomXMLFiles),
  // neither of which can change within a job, so it is built once and shared across
  // IOVs. Rebuilding it is not merely wasteful: cms::DDDetector has no destructor, so
  // dd4hep::Detector::getInstance(tag) hands back the instance built the first time and
  // the ctor's addExtension<DDVectorsMap> then re-adds an extension that is already
  // there, which segfaults inside dd4hep. That is reached whenever IdealMagneticFieldRecord
  // changes IOV, i.e. on every run transition of a data job.
  // Guarded by the kDD4hep shared resource (usesResources), which serialises these calls.
  std::shared_ptr<DDDetector> cachedXML_;
};

DDDetectorESProducer::DDDetectorESProducer(const ParameterSet& iConfig)
    : fromDB_(iConfig.getParameter<bool>("fromDB")),
      appendToDataLabel_(iConfig.getParameter<string>("appendToDataLabel")),
      confGeomXMLFiles_(fromDB_ ? "none" : iConfig.getParameter<FileInPath>("confGeomXMLFiles").fullPath()),
      rootDDName_(iConfig.getParameter<string>("rootDDName")),
      label_(iConfig.getParameter<string>("label")) {
  usesResources({edm::ESSharedResourceNames::kDD4hep});
  edm::LogVerbatim("Geometry") << "DDDetectorESProducer::fromDB " << fromDB_ << " appendToDataLabel "
                               << appendToDataLabel_ << " rootDDName " << rootDDName_ << " label " << label_
                               << " confGeomXMLFiles " << confGeomXMLFiles_;
  if (rootDDName_ == "MagneticFieldVolumes:MAGF" || rootDDName_ == "cmsMagneticField:MAGF") {
    auto c = setWhatProduced(this,
                             &DDDetectorESProducer::produceMagField,
                             edm::es::Label(iConfig.getParameter<std::string>("@module_label")));
    if (fromDB_) {
      mfToken_ = c.consumes(edm::ESInputTag("", label_));
    }
    findingRecord<IdealMagneticFieldRecord>();
  } else {
    auto c = setWhatProduced(this, &DDDetectorESProducer::produceGeom);
    if (fromDB_) {
      geomToken_ = c.consumes(edm::ESInputTag("", label_));
    }
    findingRecord<IdealGeometryRecord>();
  }
}

DDDetectorESProducer::~DDDetectorESProducer() {}

void DDDetectorESProducer::fillDescriptions(ConfigurationDescriptions& descriptions) {
  ParameterSetDescription desc;

  desc.addOptional<FileInPath>("confGeomXMLFiles");
  desc.add<string>("rootDDName", "cms:OCMS");
  desc.add<string>("label", "");
  desc.add<bool>("fromDB", false);
  descriptions.add("DDDetectorESProducer", desc);

  edm::ParameterSetDescription descDB;
  descDB.add<string>("rootDDName", "cms:OCMS");
  descDB.add<string>("label", "Extended");
  descDB.add<bool>("fromDB", true);
  descriptions.add("DDDetectorESProducerFromDB", descDB);
}

void DDDetectorESProducer::setIntervalFor(const eventsetup::EventSetupRecordKey& iKey,
                                          const IOVSyncValue& iTime,
                                          ValidityInterval& oInterval) {
  oInterval = ValidityInterval(IOVSyncValue::beginOfTime(), IOVSyncValue::endOfTime());  //infinite
}

DDDetectorESProducer::ReturnType DDDetectorESProducer::produceMagField(const IdealMagneticFieldRecord& iRecord) {
  LogVerbatim("Geometry") << "DDDetectorESProducer::Produce MF " << appendToDataLabel_;
  if (fromDB_) {
    edm::ESTransientHandle<FileBlob> gdd = iRecord.getTransientHandle(mfToken_);
    unique_ptr<vector<unsigned char> > tb = (*gdd).getUncompressedBlob();

    return make_shared<cms::DDDetector>(label_, string(tb->begin(), tb->end()), true);
  } else {
    if (!cachedXML_) {
      cachedXML_ = make_shared<DDDetector>(appendToDataLabel_, confGeomXMLFiles_);
    }
    return cachedXML_;
  }
}

DDDetectorESProducer::ReturnType DDDetectorESProducer::produceGeom(const IdealGeometryRecord& iRecord) {
  LogVerbatim("Geometry") << "DDDetectorESProducer::Produce " << appendToDataLabel_;
  if (fromDB_) {
    edm::ESTransientHandle<FileBlob> gdd = iRecord.getTransientHandle(geomToken_);
    unique_ptr<vector<unsigned char> > tb = (*gdd).getUncompressedBlob();

    return make_shared<cms::DDDetector>(label_, string(tb->begin(), tb->end()), true);
  } else {
    if (!cachedXML_) {
      cachedXML_ = make_shared<DDDetector>(appendToDataLabel_, confGeomXMLFiles_, false);
    }
    return cachedXML_;
  }
}

DDDetectorESProducer::ReturnType DDDetectorESProducer::produce() {
  LogVerbatim("Geometry") << "DDDetectorESProducer::Produce " << appendToDataLabel_;
  if (!cachedXML_) {
    cachedXML_ = make_shared<DDDetector>(appendToDataLabel_, confGeomXMLFiles_, false);
  }
  return cachedXML_;
}

DEFINE_FWK_EVENTSETUP_SOURCE(DDDetectorESProducer);
