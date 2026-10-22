// The CVH parameter catalog as a Run FlatTable, for the Runs tree of a
// calibration NanoAOD.
//
// A CVH refit maker with produceCatalog=True puts its catalog (the runtree
// columns the calibration fit reads, cvh::ParmCatalog) as an event FlatTable
// into the first event of each run that a stream processes. A stream module
// cannot put a Run product, so this module keeps the first valid table of each
// run and puts it into the Run, where the NanoAOD output writes Run FlatTables
// into the Runs tree (the label must match the keep statement, *Table). A run
// without events gets an empty table with the same columns. It must be on a
// Path: it has no event product, so in a Task (unscheduled) the framework
// would never run it on the events.

#include <memory>
#include <string>

#include "Analysis/HitAnalyzer/interface/CvhParmCatalog.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/EDPutToken.h"
#include "FWCore/Utilities/interface/InputTag.h"

class CvhParmCatalogRunTableProducer : public edm::one::EDProducer<edm::one::WatchRuns, edm::EndRunProducer> {
public:
  explicit CvhParmCatalogRunTableProducer(const edm::ParameterSet& iConfig)
      : srcToken_(consumes<nanoaod::FlatTable>(iConfig.getParameter<edm::InputTag>("src"))),
        name_(iConfig.getParameter<std::string>("name")),
        putToken_(produces<nanoaod::FlatTable, edm::Transition::EndRun>()) {}

  void beginRun(edm::Run const&, edm::EventSetup const&) override { table_.reset(); }

  void produce(edm::Event& iEvent, const edm::EventSetup&) override {
    if (table_) {
      return;
    }
    edm::Handle<nanoaod::FlatTable> handle;
    iEvent.getByToken(srcToken_, handle);
    if (handle.isValid()) {
      table_ = std::make_unique<nanoaod::FlatTable>(*handle);
    }
  }

  void endRun(edm::Run const&, edm::EventSetup const&) override {}

  void endRunProduce(edm::Run& iRun, const edm::EventSetup&) override {
    if (!table_) {
      table_ = cvh::ParmCatalog().table(name_);
    }
    iRun.put(putToken_, std::move(table_));
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag("trackrefitdimuon", "catalog"))
        ->setComment("event FlatTable of a CVH refit maker with produceCatalog=True");
    desc.add<std::string>("name", "CvhParm")->setComment("table name of the empty catalog of a run without events");
    descriptions.addWithDefaultLabel(desc);
  }

private:
  const edm::EDGetTokenT<nanoaod::FlatTable> srcToken_;
  const std::string name_;
  const edm::EDPutTokenT<nanoaod::FlatTable> putToken_;
  std::unique_ptr<nanoaod::FlatTable> table_;
};

DEFINE_FWK_MODULE(CvhParmCatalogRunTableProducer);
