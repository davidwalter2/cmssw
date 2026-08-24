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
  // Kokoulin log-T buckets for the regime-2/3 delta channel, 0 = off. The
  // default MATCHES the offline reference's own default
  // (cf_track_resolution.IONI_KOKOULIN_NBIN = 96 with IONI_KOKOULIN on), so a
  // bare run of this driver compares like with like; pass 0 to compare against
  // a python side run with IONI_KOKOULIN = 0.
  const int kokNbin = (argc > 4) ? std::atoi(argv[4]) : 96;

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

    cvhcgf::Block blk;
    blk.ioni.reserve(nstep);
    for (int i = 0; i < nstep; ++i) {
      std::istringstream rs(lines[++li]);
      // The record is 11 doubles per step historically and 13 when the export
      // ran with CVH_IONI_EXACTDELTA, which APPENDS beta2 and etot after cs so
      // that every legacy column index is unchanged. Read whatever the line
      // holds rather than a fixed count: a stride mismatch would otherwise
      // consume the next row silently.
      std::vector<double> v;
      for (double x; rs >> x;)
        v.push_back(x);
      if (v.size() < 11) {
        std::cerr << "short record row (" << v.size() << " columns)\n";
        return 1;
      }
      cvhcgf::IoniStep s;
      s.regime = static_cast<int>(v[0]);
      s.gsig2 = v[1];
      const double gam = v[9];  // Urban `scaling`, multiplies the energies
      s.a1 = v[2];
      s.e1 = v[3] * gam;
      s.a2 = v[4];
      s.e2 = v[5] * gam;
      // `a3` is a COUNT in regime 1 and xi, an ENERGY, in regime 2/3 -- so it
      // takes `scaling` in the latter and must not in the former.
      s.a3 = (s.regime >= 2) ? v[6] * gam : v[6];
      s.e0 = v[7] * gam;
      s.tmax = v[8] * gam;
      // column 10 already carries the transport weight and 1/sigma; the 1e-3
      // is the propagator's own MeV -> GeV convention on the record energies
      s.gs = v[10] * 1e-3;
      if (s.regime >= 2) {
        if (v.size() < 13) {
          std::cerr << "regime " << s.regime << " row with only " << v.size()
                    << " columns: beta2 and etot are not in the dump, and they "
                       "cannot be recovered from tmax without the particle mass\n";
          return 1;
        }
        s.beta2 = v[11];
        s.etot = v[12];
        s.kokNbin = kokNbin;
      }
      blk.ioni.push_back(s);
    }

    // Optional RADIATIVE section, written by
    // `cgf_cxx_validate.py dump --rad`:
    //
    //   RADGRID <nv> v0 v1 ... v_{nv-1}
    //   RAD <nrad>
    //   <etotMeV> <cs*w> <dEBremMeV> <dEPairMeV> <shapeBrem[nv]> <shapePair[nv]>
    //
    // The SHAPES are dumped raw and normalized here through
    // `cvhcgf::makeRadSpectrum`, so the comparison exercises the
    // normalization -- which is where the brems/pair mixture is decided --
    // and not only the transform.
    std::vector<double> vgrid;
    std::vector<std::vector<double>> radDens;
    if (li + 1 < lines.size() && lines[li + 1].rfind("RADGRID", 0) == 0) {
      std::istringstream gs(lines[++li]);
      std::string tag2;
      int nv = 0;
      gs >> tag2 >> nv;
      vgrid.resize(nv);
      for (int i = 0; i < nv; ++i)
        gs >> vgrid[i];
      std::istringstream rs(lines[++li]);
      int nrad = 0;
      rs >> tag2 >> nrad;
      radDens.reserve(nrad);
      blk.rad.reserve(nrad);
      for (int j = 0; j < nrad; ++j) {
        std::istringstream ls(lines[++li]);
        double etot = 0., csw = 0., deb = 0., dep = 0.;
        ls >> etot >> csw >> deb >> dep;
        std::vector<double> sb(nv), sp(nv);
        for (int i = 0; i < nv; ++i)
          ls >> sb[i];
        for (int i = 0; i < nv; ++i)
          ls >> sp[i];
        radDens.emplace_back(nv, 0.);
        cvhcgf::makeRadSpectrum(vgrid.data(), sb.data(), sp.data(), deb, dep, etot, nv,
                                radDens.back().data());
        cvhcgf::RadStep r;
        r.etot = etot;
        r.gs = csw * 1e-3;  // same MeV convention as IoniStep::gs
        blk.rad.push_back(r);
      }
      // the density vectors are only addressable once the vector has stopped
      // reallocating, so the pointers are attached after the loop
      for (size_t j = 0; j < blk.rad.size(); ++j)
        blk.rad[j].dNdv = radDens[j].data();
      blk.v = vgrid.data();
      blk.nv = nv;
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
      const std::complex<double> S = cvhcgf::blockExponent(blk, t);
      std::printf("SVAL %d %.17g %.17g %.17g\n", plane, t, S.real(), S.imag());
    }
    const cvhcgf::Result r = cvhcgf::inverseFisher(blk, cfg);
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
