// HadElasticTruthWatcher -- per-collision truth of the nuclear-elastic kicks a
// PRIMARY hadron takes in the tracker volume, for the survival measurement of
// the CVH nuclear-elastic family (Documents/Resolution/NUCLEAR_ELASTIC.md,
// "Survival of a kick on reconstructed tracks").
//
// WHAT IS RECORDED
// ----------------
// Fixed-length float64 records (NREC = 24 doubles), appended to `output`:
//
//   [0]  type       1 = step DEFINED by hadElastic
//                   2 = step defined by CoulombScat (single Coulomb scattering)
//                   3 = primary BEGIN (first step of the primary)
//                   4 = primary END   (the step after which the primary is not alive,
//                                      or the first step leaving the scoring volume)
//   [1]  event      G4 event id = the EDM event number
//   [2]  trackId    G4 track id (= SimTrack::trackId, PSimHit::trackId)
//   [3]  pdg
//   [4-6] x, y, z   post-step position [cm]
//   [7]  t          post-step global time [ns] -- orders the collision against the
//                   primary's PSimHits (PSimHit::tof) along the trajectory
//   [8-10] px,py,pz pre-step momentum [GeV]
//   [11] theta      angle between the pre- and post-step momentum directions [rad]
//   [12] dKE        pre - post kinetic energy [MeV]
//   [13] Tsec       kinetic energy of the secondaries THIS process created in this
//                   step [MeV] (the recoil nucleus of an elastic collision)
//   [14] secZ       atomic number of that secondary (the struck nucleus; 1 for a
//                   proton target), 0 if none
//   [15] secA       its baryon number, 0 if none
//   [16] nsec       number of secondaries THIS process created in this step
//   [17] tgtZ       Z of the nucleus the collision struck (type 1); 0 for an
//                   integral-method rejection (see realTarget)
//   [18] tgtA       its A
//   [19] steplen    [mm]
//   [20] matIdx     G4Material::GetIndex() of the pre-step material (names in
//                   `<output>.materials.txt`)
//   [21] trackLen   G4Track::GetTrackLength() after the step [cm]
//   [22] KEpre      [MeV]
//   [23] endCode    type 4: 1 hadElastic, 2 <species>Inelastic, 3 Decay,
//                   4 stopped by ionisation, 5 left the scoring volume, 6 other;
//                   types 1/2: 1 if the primary is still alive after the step
//
// DIRECTIONS (optional, `writeDirections = True`): a parallel file
// `<output>.dir` with 6 float64 per record, the pre- and post-step momentum
// DIRECTIONS (unit vectors, x y z each) -- the kick's azimuth about the
// incoming direction, which theta alone does not carry.  Row i of the
// sidecar belongs to record i of `output`.
//
// A step defined by hadElastic is NOT always a collision: G4HadronicProcess
// uses the integral method and rejects at the candidate point with probability
// 1 - xs/xs_max, returning the track unchanged (NUCLEAR_ELASTIC.md, "Geant4
// traps", 3).  `tgtZ > 0` identifies a real collision (realTarget).  The recoil
// secondary does not: a sub-threshold recoil is deposited locally.
//
// SCORING VOLUME: r < rmax, |z| < zmax (default the tracker, 120 cm x 300 cm).
// Only the primaries (ParentID == 0).  Changes no physics: the one write is the
// target-nucleus sentinel of realTarget, after the step, consuming no random
// number.
// One thread / one stream: the watcher owns its output file.

#include "SimG4Core/Watcher/interface/SimWatcher.h"
#include "SimG4Core/Watcher/interface/SimWatcherFactory.h"
#include "SimG4Core/Notification/interface/Observer.h"
#include "SimG4Core/Notification/interface/BeginOfEvent.h"
#include "SimG4Core/Notification/interface/EndOfRun.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4Event.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

class HadElasticTruthWatcher : public SimWatcher,
                               public Observer<const BeginOfEvent *>,
                               public Observer<const G4Step *>,
                               public Observer<const EndOfRun *> {
public:
  static constexpr int NREC = 24;

  explicit HadElasticTruthWatcher(const edm::ParameterSet &p)
      : out_(p.getUntrackedParameter<std::string>("output", "hadeltruth.bin")),
        rmax_(p.getUntrackedParameter<double>("rmax", 120.0)),
        zmax_(p.getUntrackedParameter<double>("zmax", 300.0)),
        writeDir_(p.getUntrackedParameter<bool>("writeDirections", false)) {
    fh_.open(out_.c_str(), std::ios::binary | std::ios::trunc);
    if (writeDir_)
      fd_.open((out_ + ".dir").c_str(), std::ios::binary | std::ios::trunc);
    std::cout << "[hadeltruth] writing " << out_ << "  rmax=" << rmax_ << " cm  zmax=" << zmax_ << " cm"
              << (fh_.is_open() ? "" : "  *** COULD NOT OPEN ***") << std::endl;
  }
  ~HadElasticTruthWatcher() override { close(); }

private:
  void update(const BeginOfEvent *evt) override {
    evtId_ = (*evt)()->GetEventID();
    begun_.clear();
    ended_.clear();
  }

  bool inside(const G4ThreeVector &x) const {
    return std::hypot(x.x(), x.y()) / CLHEP::cm < rmax_ && std::abs(x.z()) / CLHEP::cm < zmax_;
  }

  void fillCommon(const G4Step *s, int type) {
    const G4Track *t = s->GetTrack();
    const G4StepPoint *pre = s->GetPreStepPoint();
    const G4StepPoint *post = s->GetPostStepPoint();
    rec_.assign(NREC, 0.);
    rec_[0] = type;
    rec_[1] = evtId_;
    rec_[2] = t->GetTrackID();
    rec_[3] = t->GetDefinition()->GetPDGEncoding();
    const G4ThreeVector &x = post->GetPosition();
    rec_[4] = x.x() / CLHEP::cm;
    rec_[5] = x.y() / CLHEP::cm;
    rec_[6] = x.z() / CLHEP::cm;
    rec_[7] = post->GetGlobalTime() / CLHEP::ns;
    const G4ThreeVector &p0 = pre->GetMomentum();
    rec_[8] = p0.x() / CLHEP::GeV;
    rec_[9] = p0.y() / CLHEP::GeV;
    rec_[10] = p0.z() / CLHEP::GeV;
    double c = pre->GetMomentumDirection().dot(post->GetMomentumDirection());
    c = std::max(-1.0, std::min(1.0, c));
    rec_[11] = std::acos(c);
    rec_[12] = (pre->GetKineticEnergy() - post->GetKineticEnergy()) / CLHEP::MeV;
    rec_[19] = s->GetStepLength() / CLHEP::mm;
    rec_[20] = pre->GetMaterial() != nullptr ? (double)pre->GetMaterial()->GetIndex() : -1.;
    rec_[21] = t->GetTrackLength() / CLHEP::cm;
    rec_[22] = pre->GetKineticEnergy() / CLHEP::MeV;
    const G4ThreeVector &d0 = pre->GetMomentumDirection();
    const G4ThreeVector &d1 = post->GetMomentumDirection();
    dir_ = {{d0.x(), d0.y(), d0.z(), d1.x(), d1.y(), d1.z()}};
    if (pre->GetMaterial() != nullptr)
      materials_.insert(pre->GetMaterial());
  }

  void write() {
    fh_.write(reinterpret_cast<const char *>(rec_.data()), NREC * sizeof(double));
    if (writeDir_)
      fd_.write(reinterpret_cast<const char *>(dir_.data()), dir_.size() * sizeof(double));
    ++nrec_;
  }

  // The (Z, A) of the nucleus a hadElastic-defined step actually struck, or
  // (0, 1) for an integral-method rejection.  G4HadronicProcess samples the
  // target (G4CrossSectionDataStore::SampleZandA -> G4Nucleus::SetParameters)
  // only when the collision is accepted, and leaves the track unchanged
  // otherwise; the recoil secondary is no tag, since a sub-threshold recoil is
  // deposited locally.  So after every hadElastic-defined step the process's
  // target is set to the sentinel (A = 1, Z = 0) -- a free neutron, never a
  // material element -- and a step that finds the sentinel still there struck
  // nothing.  The reset happens after the step is complete and consumes no
  // random number: the physics is unchanged.
  static std::pair<int, int> realTarget(const G4Step *s) {
    const G4StepPoint *post = s->GetPostStepPoint();
    if (post == nullptr)
      return {0, 0};
    const G4VProcess *dp = post->GetProcessDefinedStep();
    if (dp == nullptr || dp->GetProcessName() != "hadElastic")
      return {0, 0};
    G4HadronicProcess *hp = const_cast<G4HadronicProcess *>(dynamic_cast<const G4HadronicProcess *>(dp));
    if (hp == nullptr)
      return {0, 0};
    G4Nucleus *n = hp->GetTargetNucleusPointer();
    const std::pair<int, int> out{n->GetZ_asInt(), n->GetA_asInt()};
    n->SetParameters(1, 0);
    return out;
  }

  void update(const G4Step *s) override {
    if (s == nullptr)
      return;
    const G4Track *t = s->GetTrack();
    if (t == nullptr)
      return;
    // the struck nucleus of THIS step (see realTarget); must run for every
    // track, since secondaries of the same species share the process instance
    const std::pair<int, int> tgt = realTarget(s);
    if (t->GetParentID() != 0)
      return;
    const int tid = t->GetTrackID();
    if (ended_.count(tid))
      return;
    const G4StepPoint *pre = s->GetPreStepPoint();
    const G4StepPoint *post = s->GetPostStepPoint();
    if (pre == nullptr || post == nullptr)
      return;

    if (!begun_.count(tid)) {
      begun_.insert(tid);
      fillCommon(s, 3);
      // BEGIN: the pre-step point
      const G4ThreeVector &x0 = pre->GetPosition();
      rec_[4] = x0.x() / CLHEP::cm;
      rec_[5] = x0.y() / CLHEP::cm;
      rec_[6] = x0.z() / CLHEP::cm;
      rec_[7] = pre->GetGlobalTime() / CLHEP::ns;
      write();
    }

    const G4VProcess *dp = post->GetProcessDefinedStep();
    const std::string pname = dp != nullptr ? std::string(dp->GetProcessName()) : std::string();
    const bool alive = t->GetTrackStatus() == fAlive;
    const bool in = inside(post->GetPosition());

    if (in && (pname == "hadElastic" || pname == "CoulombScat")) {
      fillCommon(s, pname == "hadElastic" ? 1 : 2);
      const std::vector<const G4Track *> *sec = s->GetSecondaryInCurrentStep();
      if (sec != nullptr) {
        for (const G4Track *st : *sec) {
          if (st->GetCreatorProcess() != dp)
            continue;
          rec_[16] += 1.;
          rec_[13] += st->GetKineticEnergy() / CLHEP::MeV;
          const G4ParticleDefinition *d = st->GetDefinition();
          if (d->GetBaryonNumber() > rec_[15]) {
            rec_[14] = d->GetAtomicNumber() > 0 ? d->GetAtomicNumber() : (d->GetPDGEncoding() == 2212 ? 1 : 0);
            rec_[15] = d->GetBaryonNumber();
          }
        }
      }
      if (pname == "hadElastic") {
        rec_[17] = tgt.first;
        rec_[18] = tgt.second;
      }
      rec_[23] = alive ? 1. : 0.;
      write();
    }

    if (!alive || !in) {
      ended_.insert(tid);
      fillCommon(s, 4);
      int code = 6;
      if (!in)
        code = 5;
      else if (pname == "hadElastic")
        code = 1;
      else if (pname.find("Inelastic") != std::string::npos)
        code = 2;
      else if (pname == "Decay")
        code = 3;
      else if (post->GetKineticEnergy() <= 0.)
        code = 4;
      rec_[23] = code;
      write();
    }
  }

  void update(const EndOfRun *) override { close(); }

  void close() {
    if (closed_)
      return;
    closed_ = true;
    if (fh_.is_open())
      fh_.close();
    if (fd_.is_open())
      fd_.close();
    std::ofstream mf((out_ + ".materials.txt").c_str(), std::ios::trunc);
    for (const G4Material *m : materials_)
      mf << m->GetIndex() << " " << m->GetName() << " " << m->GetDensity() / (CLHEP::g / CLHEP::cm3) << "\n";
    std::cout << "[hadeltruth] wrote " << nrec_ << " records of " << NREC << " float64 to " << out_ << std::endl;
  }

  const std::string out_;
  const double rmax_, zmax_;
  std::ofstream fh_;
  const bool writeDir_;
  std::ofstream fd_;
  std::array<double, 6> dir_{};
  bool closed_ = false;
  long nrec_ = 0;
  double evtId_ = -1;
  std::set<int> begun_, ended_;
  std::set<const G4Material *> materials_;
  std::vector<double> rec_ = std::vector<double>(NREC, 0.);
};

DEFINE_SIMWATCHER(HadElasticTruthWatcher);
