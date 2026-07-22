#include "TrackPropagation/Geant4e/interface/MaterialGroupModel.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "G4LogicalVolume.hh"

#include <fstream>
#include <sstream>

MaterialGroupModel::MaterialGroupModel(const std::string &rulesFile) {
  std::ifstream in(rulesFile);
  if (!in) {
    throw cms::Exception("MaterialGroupModel") << "cannot open rules file " << rulesFile;
  }
  // group 0 = catch-all
  gnames_ = {"other"};
  kval_ = {0.};
  prior_ = {0.2};

  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') {
      continue;
    }
    std::istringstream ss(line);
    std::string tag, gname, rex, srmin, srmax, szmin, szmax;
    int gid, zside;
    double kinit, prior;
    ss >> tag;
    if (tag != "RULE") {
      continue;
    }
    if (!(ss >> gid >> gname >> rex >> srmin >> srmax >> szmin >> szmax >> zside >> kinit >> prior)) {
      throw cms::Exception("MaterialGroupModel") << "malformed RULE line: " << line;
    }
    auto bound = [](const std::string &s) { return s == "-" ? -1. : std::stod(s); };
    Rule r;
    r.groupId = gid;
    r.rex = std::regex(rex);
    r.rmin = bound(srmin);
    r.rmax = bound(srmax);
    r.zmin = bound(szmin);
    r.zmax = bound(szmax);
    r.zside = zside;
    rules_.push_back(std::move(r));

    if (gid >= static_cast<int>(gnames_.size())) {
      gnames_.resize(gid + 1);
      kval_.resize(gid + 1, 0.);
      prior_.resize(gid + 1, 0.2);
    }
    gnames_[gid] = gname;
    kval_[gid] = kinit;
    prior_[gid] = prior;
  }
  nGroups_ = gnames_.size();
  if (rules_.empty()) {
    throw cms::Exception("MaterialGroupModel") << "no RULE lines in " << rulesFile;
  }
}

const std::vector<unsigned short> &MaterialGroupModel::shortlist(const G4LogicalVolume *lv) const {
  auto it = lvCache_.find(lv);
  if (it != lvCache_.end()) {
    return it->second;
  }
  std::vector<unsigned short> sl;
  const std::string name = lv->GetName();
  for (unsigned short i = 0; i < rules_.size(); ++i) {
    if (std::regex_search(name, rules_[i].rex)) {
      sl.push_back(i);
    }
  }
  return lvCache_.emplace(lv, std::move(sl)).first->second;
}

int MaterialGroupModel::classify(const G4LogicalVolume *lv, double r_cm, double z_cm) const {
  for (unsigned short i : shortlist(lv)) {
    const Rule &r = rules_[i];
    if (r.rmin >= 0. && r_cm < r.rmin) {
      continue;
    }
    if (r.rmax >= 0. && r_cm >= r.rmax) {
      continue;
    }
    const double az = std::abs(z_cm);
    if (r.zmin >= 0. && az < r.zmin) {
      continue;
    }
    if (r.zmax >= 0. && az >= r.zmax) {
      continue;
    }
    if (r.zside != 0 && (z_cm < 0. ? -1 : +1) != r.zside) {
      continue;
    }
    return r.groupId;
  }
  return 0;
}

double MaterialGroupModel::materialOffset(const G4LogicalVolume *lv, double r_cm, double z_cm) const {
  const int g = classify(lv, r_cm, z_cm);
  double k = kval_[g];
  if (g == injGroup_) {
    k += injEps_;
  }
  return k;
}
