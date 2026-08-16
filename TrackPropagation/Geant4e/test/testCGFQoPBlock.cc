// Standalone check of cvhcgf against the validated python reference.
//
// Reads the dump written by
//   calibration_studies/resolution/cgf_cxx_validate.py dump
// and prints S(t) and 1/I in a form that script's `cmp` mode reads back.
//
// Deliberately free of every CMSSW and Geant4 dependency so it can be built
// with a bare g++ and so a failure here is unambiguously the CGF code's:
//   g++ -O2 -std=c++17 -I../interface testCGFQoPBlock.cc ../src/CGFQoPBlock.cc -o testcgf

#include "TrackPropagation/Geant4e/interface/CGFQoPBlock.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "usage: testcgf <dumpfile> [nt] [npad]\n";
    return 1;
  }
  cvhcgf::Config cfg;
  if (argc > 2)
    cfg.nt = std::atoi(argv[2]);
  if (argc > 3)
    cfg.npad = std::atoi(argv[3]);

  std::ifstream in(argv[1]);
  if (!in) {
    std::cerr << "cannot open " << argv[1] << "\n";
    return 1;
  }

  std::vector<std::string> lines;
  {
    std::string l;
    while (std::getline(in, l))
      lines.push_back(l);
  }

  for (size_t li = 0; li < lines.size(); ++li) {
    const std::string &line = lines[li];
    if (line.rfind("BLOCK", 0) != 0)
      continue;
    std::istringstream hs(line);
    std::string tag;
    int plane = 0, nstep = 0;
    double sigma = 0., invIref = 0., topref = 0.;
    hs >> tag >> plane >> nstep >> sigma >> invIref >> topref;

    std::vector<cvhcgf::IoniStep> steps;
    steps.reserve(nstep);
    for (int i = 0; i < nstep; ++i) {
      std::istringstream rs(lines[++li]);
      double v[11];
      for (double &x : v)
        rs >> x;
      cvhcgf::IoniStep s;
      s.regime = static_cast<int>(v[0]);
      s.gsig2 = v[1];
      const double gam = v[9];  // Urban `scaling`, multiplies the energies
      s.a1 = v[2];
      s.e1 = v[3] * gam;
      s.a2 = v[4];
      s.e2 = v[5] * gam;
      s.a3 = v[6];
      s.e0 = v[7] * gam;
      s.tmax = v[8] * gam;
      // column 10 already carries the transport weight and 1/sigma; the 1e-3
      // is the propagator's own MeV -> GeV convention on the record energies
      s.gs = v[10] * 1e-3;
      steps.push_back(s);
    }

    ++li;  // SREF
    std::vector<double> tl;
    while (li + 1 < lines.size() && lines[li + 1].rfind("BLOCK", 0) != 0 && !lines[li + 1].empty()) {
      std::istringstream ss(lines[++li]);
      double t, a, b;
      ss >> t >> a >> b;
      tl.push_back(t);
    }

    for (double t : tl) {
      const std::complex<double> S = cvhcgf::blockExponent(steps, t);
      std::printf("SVAL %d %.17g %.17g %.17g\n", plane, t, S.real(), S.imag());
    }
    const cvhcgf::Result r = cvhcgf::inverseFisher(steps, cfg);
    std::printf("IVAL %d %.17g\n", plane, r.invFisher);
    std::printf("MVAL %d %.17g %d\n", plane, r.zmode, r.nZeroCross);
    // score at a fixed ladder of residuals, spanning the core and both tails
    for (double z : {-20., -10., -5., -3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 8., 10.}) {
      bool cl = false;
      const double ps = cvhcgf::scoreAt(r, z, &cl);
      std::printf("PVAL %d %.17g %.17g %d\n", plane, z, ps, int(cl));
    }
    std::printf("# plane %d ok=%d nstep=%d kappa2=%.8g tmax=%.6f mass=%.8f "
                "massfrac=%.8f zmode=%.5f support=[%.4f,%.4f] pyref=%.8f\n",
                plane, int(r.ok), r.nsteps, r.kappa2, r.tmax, r.mass, r.massFrac, r.zmode, r.zlo, r.zhi,
                invIref);
  }
  return 0;
}
