#ifndef Analysis_HitAnalyzer_CvhParmCatalog_h
#define Analysis_HitAnalyzer_CvhParmCatalog_h

// The catalog of the global correction parameters of a CVH refit: one row per
// parameter, in global-index order, with the columns of the `runtree` that the
// calibration fit reads (parmtype, rawdetid, subdet, layer for the priors and
// the parameter blocks; position and eta/phi to locate a module). As a
// NanoAOD FlatTable it lands in the Runs tree, so a calibration NanoAOD
// carries the meaning of the indices its gradients and Hessians refer to.
// The column set is defined here once, so that an empty catalog (a run without
// events) has the same columns as a filled one.

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "DataFormats/NanoAOD/interface/FlatTable.h"

namespace cvh {

  struct ParmCatalog {
    std::vector<int> parmtype, subdet, layer, stereo, glued;
    std::vector<uint32_t> rawdetid;
    std::vector<float> x, y, z, eta, phi;

    void clear() {
      parmtype.clear();
      subdet.clear();
      layer.clear();
      stereo.clear();
      glued.clear();
      rawdetid.clear();
      x.clear();
      y.clear();
      z.clear();
      eta.clear();
      phi.clear();
    }

    void reserve(std::size_t n) {
      parmtype.reserve(n);
      subdet.reserve(n);
      layer.reserve(n);
      stereo.reserve(n);
      glued.reserve(n);
      rawdetid.reserve(n);
      x.reserve(n);
      y.reserve(n);
      z.reserve(n);
      eta.reserve(n);
      phi.reserve(n);
    }

    void push_back(int iparmtype,
                   uint32_t irawdetid,
                   int isubdet,
                   int ilayer,
                   int istereo,
                   int iglued,
                   double ix,
                   double iy,
                   double iz,
                   double ieta,
                   double iphi) {
      parmtype.push_back(iparmtype);
      rawdetid.push_back(irawdetid);
      subdet.push_back(isubdet);
      layer.push_back(ilayer);
      stereo.push_back(istereo);
      glued.push_back(iglued);
      x.push_back(ix);
      y.push_back(iy);
      z.push_back(iz);
      eta.push_back(ieta);
      phi.push_back(iphi);
    }

    std::size_t size() const { return parmtype.size(); }

    std::unique_ptr<nanoaod::FlatTable> table(const std::string& name) const {
      auto t = std::make_unique<nanoaod::FlatTable>(size(), name, false, false);
      t->setDoc("CVH global correction parameters, one row per global index (the runtree of the refit)");
      t->addColumn<int>("parmtype", parmtype, "parameter type");
      t->addColumn<uint32_t>("rawdetid", rawdetid, "raw DetId of the module; mode or group index of a global block");
      t->addColumn<int>("subdet", subdet, "tracker subdetector (-99 for global blocks)");
      t->addColumn<int>("layer", layer, "layer (signed, -99 for global blocks)");
      t->addColumn<int>("stereo", stereo, "stereo module");
      t->addColumn<int>("glued", glued, "glued module");
      t->addColumn<float>("x", x, "module position x (cm)", -1);
      t->addColumn<float>("y", y, "module position y (cm)", -1);
      t->addColumn<float>("z", z, "module position z (cm)", -1);
      t->addColumn<float>("eta", eta, "module eta", -1);
      t->addColumn<float>("phi", phi, "module phi", -1);
      return t;
    }
  };

}  // namespace cvh

#endif
