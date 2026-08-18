// Deactivates named Geant4 processes on the TRACKING thread, and counts how
// many steps each process actually defined, so that "the switch took" is a
// measurement rather than an assumption.
//
// THREE ROUTES THAT LOOK RIGHT.  TWO OF THEM ARE SILENTLY WRONG.
// --------------------------------------------------------------
// 1. `process.g4SimHits.G4Commands = ['/process/inactivate muBrems']`
//    INERT.  CMSSW applies G4Commands from RunManagerMT::initG4 and
//    RunManagerMTWorker::initializeG4 while Geant4 is still in G4State_PreInit
//    (both call G4UImanager::ApplyCommand and only then SetNewState(Init)),
//    i.e. before InitializePhysics has constructed any process.  ApplyCommand
//    returns an error and CMSSW DISCARDS the return code: nothing printed,
//    nothing done.  Measured: two 2000-event pT = 40 runs at a fixed seed, one
//    with those commands and one without, came out BIT-IDENTICAL event for
//    event (0 of 2000 differing).
//
// 2. A watcher living OUTSIDE the Simulation biglib (e.g. in
//    Analysis/HitAnalyzer/plugins), acting through
//    G4ParticleTable::FindParticle("mu-") or ...->GetProcessManager().
//    ALSO INERT, and worse.  src/BigProducts/Simulation/BuildFile.xml carries
//    <use name="geant4static"/> + DROP_DEP="geant4core", so pluginSimulation.so
//    -- where OscarMTProducer and the whole G4 machinery live -- contains a
//    PRIVATE, statically linked Geant4 with local symbols
//    (nm -C --defined-only shows 't G4ProcessTable::GetProcessTable()').  A
//    plugin in lib/ links the SHARED libG4*.so and therefore talks to a second,
//    never-initialised Geant4: FindParticle returns nullptr,
//    G4ProcessTable::Length() is 0, and GetProcessManager() SEGFAULTS on the
//    tracking thread.  Objects passed in by the simulation (G4Step, G4Track)
//    are fine; every G4 SINGLETON is the wrong one.  This file therefore has to
//    live in a package that is part of BigProducts/Simulation.
//
// 3. What is used here: G4ProcessTable::SetProcessActivation(name, false) from
//    BeginOfTrack.  G4ProcessTable is thread-local and its elements store the
//    G4ProcessManager pointer that was registered WITH the process, so it never
//    has to resolve a manager from a particle definition -- which is the call
//    that segfaults.  BeginOfTrack runs on the tracking thread and before the
//    track is stepped.
//
// THE EVIDENCE IS THE STEP COUNTER, NOT THE FLAG
// ----------------------------------------------
// Every step is attributed to the process that defined it
// (GetPostStepPoint()->GetProcessDefinedStep()).  At EndOfRun the per-process
// counts are printed for the primary and for all tracks.  A working
// deactivation shows EXACTLY ZERO muBrems and muPairProd steps while muIoni is
// unchanged; a failed one shows the ordinary non-zero counts.  That is a direct
// statement about what Geant4 ran, independent of any flag.
//
// SCOPE: additive and default-inert.  With an empty `inactivate` list this only
// counts and prints.

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4EmParameters.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronicInteraction.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessTable.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4Threading.hh"
#include "G4ios.hh"
#include <atomic>
#include <cstdlib>
#include <iostream>

#include <map>
#include <string>
#include <vector>

class ProcessActivationWatcher : public SimWatcher,
                                 public Observer<const BeginOfTrack *>,
                                 public Observer<const G4Step *>,
                                 public Observer<const EndOfRun *> {
public:
  explicit ProcessActivationWatcher(const edm::ParameterSet &p);
  ~ProcessActivationWatcher() override { report("dtor"); }

private:
  void update(const BeginOfTrack *) override;
  void update(const G4Step *) override;
  void update(const EndOfRun *) override;
  void report(const char *where);

  std::vector<std::string> inactivate_;
  std::vector<std::string> activate_;
  const bool dumpHadModels_;
  const bool dumpEmParams_;
  bool applied_ = false;
  bool reported_ = false;
  std::map<std::string, long> nprim_, nall_;
};

ProcessActivationWatcher::ProcessActivationWatcher(const edm::ParameterSet &p)
    : inactivate_(
          p.getUntrackedParameter<std::vector<std::string>>("inactivate", std::vector<std::string>())),
      activate_(p.getUntrackedParameter<std::vector<std::string>>("activate", std::vector<std::string>())),
      // ParameterSet, not getenv: this watcher already had a PSet, and a
      // diagnostic that leaves no trace in the job configuration cannot be
      // matched to the output it produced.
      dumpHadModels_(p.getUntrackedParameter<bool>("dumpHadronicModels", false)),
      dumpEmParams_(p.getUntrackedParameter<bool>("dumpEmParameters", false)) {
  std::cout << "[procact] constructed: inactivate=" << inactivate_.size() << " activate=" << activate_.size()
         << std::endl;
  for (const auto &n : inactivate_)
    std::cout << "[procact]   requested INACTIVE '" << n << "'" << std::endl;
  for (const auto &n : activate_)
    std::cout << "[procact]   requested ACTIVE   '" << n << "'" << std::endl;
}

void ProcessActivationWatcher::update(const BeginOfTrack *) {
  // dumpHadronicModels: which hadronic model and cross-section set actually
  // handles each hadronic process for THIS particle, in THIS physics list.
  // "Following the G4 implementation" is only meaningful against the concrete
  // assignment, and FTFP_BERT_EMM picks different models per species and
  // energy range; reading it off the process table is authoritative where
  // reading the physics-list source is not (CMS's source is not on cvmfs).
  {
    static std::atomic<bool> dumpedHad{false};
    bool exp2 = false;
    if (dumpHadModels_ && dumpedHad.compare_exchange_strong(exp2, true)) {
      std::cout << "### CVH_HADMODELS_BEGIN" << std::endl;
      // per PARTICLE, which is what determines the assignment we must follow
      auto *ptab = G4ParticleTable::GetParticleTable();
      for (G4int i = 0; i < ptab->size(); ++i) {
        G4ParticleDefinition *pd = ptab->GetParticle(i);
        if (pd == nullptr) continue;
        const G4String pn = pd->GetParticleName();
        if (pn != "pi+" && pn != "pi-" && pn != "kaon+" && pn != "kaon-" &&
            pn != "proton" && pn != "anti_proton" && pn != "mu+" && pn != "mu-")
          continue;
        G4ProcessManager *pm = pd->GetProcessManager();
        if (pm == nullptr) continue;
        G4ProcessVector *pv = pm->GetProcessList();
        for (G4int k = 0; k < (G4int)pv->size(); ++k) {
          auto *hp = dynamic_cast<G4HadronicProcess *>((*pv)[k]);
          if (hp == nullptr) continue;
          for (auto *hi : hp->GetHadronicInteractionList()) {
            if (hi == nullptr) continue;
            std::cout << "  " << pn << "  " << hp->GetProcessName() << "  -> " << hi->GetModelName()
                      << "  E=[" << hi->GetMinEnergy()/CLHEP::GeV << ","
                      << hi->GetMaxEnergy()/CLHEP::GeV << "] GeV" << std::endl;
          }
        }
      }
      std::cout << "### CVH_HADMODELS_END" << std::endl;
    }
  }

  // dumpEmParameters: the EM parameter block, once. G4EmParameters is a
  // GLOBAL SINGLETON, and the sim and the model are separate jobs running
  // different physics lists -- the sim CMS's, the model
  // G4ErrorPhysicsListForCVH, which sets no EM parameter at all. So anything
  // the model reads from this singleton at Initialise() time silently takes a
  // Geant4 default where the sim takes a CMS value. That has bitten here
  // before (a fork queried Spline() and got a different default, see
  // G4TablesForExtrapolatorForCVH.h). Dumped here rather than at BeginOfRun
  // because the EM tables are built lazily and are not complete until the
  // first track.
  {
    static std::atomic<bool> dumped{false};
    bool expected = false;
    if (dumpEmParams_ && dumped.compare_exchange_strong(expected, true)) {
      // std::cout, NOT G4cout: SimG4Core installs a G4UIsession that captures
      // G4cout, so a G4cout dump here is swallowed while the same dump in the
      // model job (no such session) appears. That difference cost a run.
      std::cout << "### CVH_EMPARAMS_BEGIN tag=SIM" << std::endl;
      G4EmParameters::Instance()->StreamInfo(std::cout);
      std::cout << "### CVH_EMPARAMS_END" << std::endl;
    }
  }

  if (applied_)
    return;
  applied_ = true;
  G4ProcessTable *tbl = G4ProcessTable::GetProcessTable();
  std::cout << "[procact] BeginOfTrack tid=" << G4Threading::G4GetThreadId()
         << ": process table length=" << (tbl != nullptr ? (int)tbl->Length() : -1) << std::endl;
  if (tbl == nullptr)
    return;
  for (const auto &n : inactivate_) {
    tbl->SetProcessActivation(G4String(n), false);
    std::cout << "[procact]   SetProcessActivation(" << n << ", false)" << std::endl;
  }
  for (const auto &n : activate_) {
    tbl->SetProcessActivation(G4String(n), true);
    std::cout << "[procact]   SetProcessActivation(" << n << ", true)" << std::endl;
  }
}

void ProcessActivationWatcher::update(const G4Step *step) {
  const G4VProcess *pr = step->GetPostStepPoint()->GetProcessDefinedStep();
  const std::string nm = (pr != nullptr) ? std::string(pr->GetProcessName()) : std::string("none");
  ++nall_[nm];
  if (step->GetTrack()->GetParentID() == 0)
    ++nprim_[nm];
}

void ProcessActivationWatcher::update(const EndOfRun *) { report("EndOfRun"); }

void ProcessActivationWatcher::report(const char *where) {
  if (reported_)
    return;
  reported_ = true;
  std::cout << "[procact] " << where << ": steps by defining process "
         << "(a deactivated process must show EXACTLY 0)" << std::endl;
  for (const auto &kv : nall_)
    std::cout << "[procact]   " << kv.first << "  primary " << nprim_[kv.first] << "  all " << kv.second
           << std::endl;
}

DEFINE_SIMWATCHER(ProcessActivationWatcher);
