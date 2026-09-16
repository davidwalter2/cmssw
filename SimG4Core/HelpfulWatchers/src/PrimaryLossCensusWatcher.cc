// Per-event census of WHICH GEANT4 PROCESS TAKES THE PRIMARY'S ENERGY, so that
// "something lives in the tail of the q/p residual" can be converted into a
// named process.
//
// WHY IT IS HERE AND NOT IN Analysis/HitAnalyzer/plugins
// ------------------------------------------------------
// It dumps G4ProductionCutsTable, which is a G4 SINGLETON.
// src/BigProducts/Simulation/BuildFile.xml carries <use name="geant4static"/>
// + DROP_DEP="geant4core", so pluginSimulation.so owns a PRIVATE, statically
// linked Geant4 with local symbols; a plugin living in lib/ links the shared
// libG4*.so and therefore talks to a second, never-initialised Geant4 (empty
// particle table, zero-length process table, segfault on
// GetProcessManager()).  SimG4Core/HelpfulWatchers IS in the
// BigProducts/Simulation list, so a plugin added here is compiled INTO
// pluginSimulation.so and sees the simulation's own Geant4.  This is the same
// constraint ProcessActivationWatcher.cc documents.
//
// WHY A RAW BINARY FILE AND NOT ROOT
// ----------------------------------
// SimG4Core/HelpfulWatchers/BuildFile.xml has no ROOT dependency and adding one
// relinks the whole simulation biglib for no physics reason.  The record is a
// fixed-length float32 array per event, written with std::ofstream; the offline
// side reads it with numpy in one np.fromfile.  ONE RECORD PER EVENT, written
// at EndOfEvent unconditionally, in the same order ToyStateNtuplizer fills its
// tree -- so entry i here is entry i there and the join is by index.  Events in
// which the primary never enters the scoring region still get a record (all
// zeros but ekin0), which is what keeps the two files aligned.
//
// WHAT IS COUNTED
// ---------------
// Only the PRIMARY (ParentID == 0), and only while its pre-step radius is below
// `rmax` -- the outermost scoring shell -- so the census covers exactly the
// path the closure is evaluated over and not the calorimeter beyond it.
//
// For every step:  dKE = preKE - postKE  is attributed to the process that
// DEFINED the step (GetPostStepPoint()->GetProcessDefinedStep()).  A step
// defined by Transportation or by a boundary still carries the continuous
// along-step ionization loss, which is why the per-process sums are of dKE and
// not of the discrete part alone.  Separately, the energy handed to explicit
// SECONDARIES in that step (GetSecondaryInCurrentStep()) is summed per process:
// that is the delta-ray / brems / pair energy that Geant4 puts on a new track
// instead of into the continuous straggling integral.  dKE - Esec is the
// continuous part.  The split between the two is exactly what a production-cut
// mismatch moves, so both are recorded.
//
// SCOPE: read-only.  It changes no physics and, with no parameters set, only
// writes a file.

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/BeginOfTrack.h"
#include "SimG4Core/Notification/interface/EndOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4ProductionCutsTable.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {
  // FIXED process codes.  A name->id map built at run time would make the
  // offline reader depend on parsing the log; this cannot go stale.
  enum ProcCode {
    kNone = 0,
    kMuIoni = 1,
    kMuBrems = 2,
    kMuPairProd = 3,
    kMuonNuclear = 4,
    kCoulombScat = 5,
    kMsc = 6,
    kTransportation = 7,
    kDecay = 8,
    kOtherIoni = 9,
    kOther = 10,
    kNProc = 11
  };

  int procCode(const G4VProcess *p) {
    if (p == nullptr)
      return kNone;
    const G4String &n = p->GetProcessName();
    if (n == "muIoni")
      return kMuIoni;
    if (n == "muBrems")
      return kMuBrems;
    if (n == "muPairProd")
      return kMuPairProd;
    if (n == "muonNuclear")
      return kMuonNuclear;
    if (n == "CoulombScat")
      return kCoulombScat;
    if (n == "msc" || n == "TransportationWithMsc")
      return kMsc;
    if (n == "Transportation")
      return kTransportation;
    if (n == "Decay" || n == "DecayWithSpin")
      return kDecay;
    if (n == "hIoni" || n == "eIoni" || n == "ionIoni")
      return kOtherIoni;
    return kOther;
  }

  // record layout, float32, NREC entries per event -- mirrored in the offline
  // reader (calibration_studies/resolution/tail_probe.py: REC_FIELDS)
  enum {
    iEkin0 = 0,
    iNstep = 1,
    iDEtot = 2,
    iDEdep = 3,
    iDEsec = 4,
    iDEproc = 5,               // .. 5 + kNProc - 1
    iEsecProc = iDEproc + kNProc,   // .. + kNProc - 1
    iNsecProc = iEsecProc + kNProc,
    iMaxdE = iNsecProc + kNProc,
    iMaxdECode,
    iMaxdER,
    iMaxSec,
    iMaxSecCode,
    iMaxSecR,
    iNhard,
    iRlast,
    NREC
  };
}  // namespace

class PrimaryLossCensusWatcher : public SimWatcher,
                                 public Observer<const BeginOfTrack *>,
                                 public Observer<const BeginOfEvent *>,
                                 public Observer<const G4Step *>,
                                 public Observer<const EndOfEvent *>,
                                 public Observer<const EndOfRun *> {
public:
  explicit PrimaryLossCensusWatcher(const edm::ParameterSet &p);
  ~PrimaryLossCensusWatcher() override { closeOut(); }

private:
  void update(const BeginOfTrack *) override;
  void update(const BeginOfEvent *) override;
  void update(const G4Step *) override;
  void update(const EndOfEvent *) override;
  void update(const EndOfRun *) override;
  void dumpCuts();
  void closeOut();

  double rmax_;      // cm
  double hard_;      // MeV, "hard step" threshold
  std::string out_;
  std::ofstream fh_;
  bool cutsDumped_ = false;
  bool closed_ = false;
  long nev_ = 0;
  std::vector<float> rec_;
};

PrimaryLossCensusWatcher::PrimaryLossCensusWatcher(const edm::ParameterSet &p)
    : rmax_(p.getUntrackedParameter<double>("rmax", 107.0)),
      hard_(p.getUntrackedParameter<double>("hard", 1.0)),
      out_(p.getUntrackedParameter<std::string>("output", "primaryloss.bin")),
      rec_(NREC, 0.f) {
  fh_.open(out_.c_str(), std::ios::binary | std::ios::trunc);
  std::cout << "[plcensus] rmax=" << rmax_ << " cm  hard=" << hard_ << " MeV  NREC=" << NREC << "  -> " << out_
            << (fh_.is_open() ? "  (open)" : "  *** COULD NOT OPEN ***") << std::endl;
}

void PrimaryLossCensusWatcher::update(const BeginOfTrack *) {
  if (!cutsDumped_) {
    cutsDumped_ = true;
    dumpCuts();
  }
}

void PrimaryLossCensusWatcher::dumpCuts() {
  G4ProductionCutsTable *t = G4ProductionCutsTable::GetProductionCutsTable();
  if (t == nullptr) {
    std::cout << "[plcensus] no G4ProductionCutsTable" << std::endl;
    return;
  }
  const std::size_t n = t->GetTableSize();
  std::cout << "[plcensus] G4ProductionCutsTable: " << n << " couples.  "
            << "Delta rays BELOW the e- energy cut are NOT made as secondaries: their "
               "energy is lost continuously, i.e. the SIM's straggling integral is "
               "RESTRICTED to tcut = that energy, while the offline record is exported "
               "with tcut = Tmax."
            << std::endl;
  std::cout << "[plcensus]   idx  material                 rho[g/cm3]   "
               "range[mm]    Ecut_gamma[MeV]  Ecut_e-[MeV]   Ecut_e+[MeV]"
            << std::endl;
  for (std::size_t i = 0; i < n; ++i) {
    const G4MaterialCutsCouple *c = t->GetMaterialCutsCouple((G4int)i);
    if (c == nullptr || c->GetMaterial() == nullptr)
      continue;
    const G4Material *m = c->GetMaterial();
    const std::vector<G4double> *eg = t->GetEnergyCutsVector(0);
    const std::vector<G4double> *ee = t->GetEnergyCutsVector(1);
    const std::vector<G4double> *ep = t->GetEnergyCutsVector(2);
    const G4ProductionCuts *pc = c->GetProductionCuts();
    std::cout << "[plcensus]   " << i << "  " << m->GetName() << "   " << m->GetDensity() / (g / cm3) << "   "
              << (pc != nullptr ? pc->GetProductionCut(1) / mm : -1.) << "   "
              << (eg != nullptr && i < eg->size() ? (*eg)[i] / MeV : -1.) << "   "
              << (ee != nullptr && i < ee->size() ? (*ee)[i] / MeV : -1.) << "   "
              << (ep != nullptr && i < ep->size() ? (*ep)[i] / MeV : -1.) << std::endl;
  }
}

void PrimaryLossCensusWatcher::update(const BeginOfEvent *) { rec_.assign(NREC, 0.f); }

void PrimaryLossCensusWatcher::update(const G4Step *step) {
  const G4Track *trk = step->GetTrack();
  if (trk->GetParentID() != 0)
    return;
  const G4StepPoint *pre = step->GetPreStepPoint();
  const G4StepPoint *post = step->GetPostStepPoint();
  const G4ThreeVector pp = pre->GetPosition() / CLHEP::cm;
  const double r = std::hypot(pp.x(), pp.y());
  if (r > rmax_)
    return;

  if (rec_[iEkin0] == 0.f)
    rec_[iEkin0] = (float)(pre->GetKineticEnergy() / CLHEP::MeV);

  const double dKE = (pre->GetKineticEnergy() - post->GetKineticEnergy()) / CLHEP::MeV;
  const int code = procCode(post->GetProcessDefinedStep());

  double esec = 0.;
  double emax = 0.;
  const std::vector<const G4Track *> *sec = step->GetSecondaryInCurrentStep();
  if (sec != nullptr) {
    for (const G4Track *s : *sec) {
      const double e = s->GetKineticEnergy() / CLHEP::MeV;
      esec += e;
      if (e > emax)
        emax = e;
      rec_[iNsecProc + code] += 1.f;
    }
  }

  rec_[iNstep] += 1.f;
  rec_[iDEtot] += (float)dKE;
  rec_[iDEdep] += (float)(step->GetTotalEnergyDeposit() / CLHEP::MeV);
  rec_[iDEsec] += (float)esec;
  rec_[iDEproc + code] += (float)dKE;
  rec_[iEsecProc + code] += (float)esec;
  rec_[iRlast] = (float)r;

  if (dKE > rec_[iMaxdE]) {
    rec_[iMaxdE] = (float)dKE;
    rec_[iMaxdECode] = (float)code;
    rec_[iMaxdER] = (float)r;
  }
  if (emax > rec_[iMaxSec]) {
    rec_[iMaxSec] = (float)emax;
    rec_[iMaxSecCode] = (float)code;
    rec_[iMaxSecR] = (float)r;
  }
  if (dKE > hard_)
    rec_[iNhard] += 1.f;
}

void PrimaryLossCensusWatcher::update(const EndOfEvent *) {
  if (fh_.is_open())
    fh_.write(reinterpret_cast<const char *>(rec_.data()), (std::streamsize)(NREC * sizeof(float)));
  ++nev_;
}

void PrimaryLossCensusWatcher::update(const EndOfRun *) { closeOut(); }

void PrimaryLossCensusWatcher::closeOut() {
  if (closed_)
    return;
  closed_ = true;
  if (fh_.is_open()) {
    fh_.flush();
    fh_.close();
  }
  std::cout << "[plcensus] wrote " << nev_ << " records of " << NREC << " float32 to " << out_ << std::endl;
}

DEFINE_SIMWATCHER(PrimaryLossCensusWatcher);
