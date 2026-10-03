// Standalone smoke test for ScalarPot3DEval.
//
// Reads the dump file path from argv[1], evaluates basis values and
// the absolute field at a small set of test points, and prints them
// in a stable text format consumed by the mfs-side comparison driver
// test_scalar_pot3d_against_mfs.py.
//
// Build: scram b in MagneticField/ParametrizedEngine, then run via
//   testScalarPot3DEval <dump.txt> > cpp_out.txt
//
// Run inside the el7 container (see CLAUDE.md).

#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "MagneticField/ParametrizedEngine/interface/ScalarPot3DEval.h"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <dump.txt>\n", argv[0]);
    return 1;
  }
  magfieldparam::ScalarPot3DEval ev(argv[1]);

  std::mt19937 rng(20260508);
  std::uniform_real_distribution<double> ur(5.0, 290.0);
  std::uniform_real_distribution<double> uphi(-3.14159265358979, 3.14159265358979);
  std::uniform_real_distribution<double> uz(-300.0, 300.0);

  const unsigned int N_PTS = 64;
  std::vector<double> bz, br, bphi, bx, by;

  std::cout << std::scientific << std::setprecision(16);
  // Header so the Python driver can sanity-check what was tested
  std::cout << "# nModes " << ev.nModes() << "\n";
  std::cout << "# r_scale " << ev.rScale() << "\n";
  std::cout << "# z0 " << ev.z0() << "\n";
  std::cout << "# cmssw_norm " << (ev.cmsswNorm() ? 1 : 0) << "\n";
  std::cout << "# l_max " << ev.lMax() << "\n";
  std::cout << "# n_test_points " << N_PTS << "\n";
  std::cout << "# columns: x_cm y_cm z_cm  Bx_T By_T Bz_T  basis_Bz_first basis_Br_first basis_Bphi_first basis_Bz_last basis_Br_last basis_Bphi_last\n";

  for (unsigned int k = 0; k < N_PTS; ++k) {
    const double r = ur(rng);
    const double phi = uphi(rng);
    const double z = uz(rng);
    const double x = r * std::cos(phi);
    const double y = r * std::sin(phi);
    const GlobalPoint gp(x, y, z);

    const GlobalVector B = ev.evaluateAbsoluteAt(gp, ev.initCoeffs());
    ev.evaluateBasisAt(gp, bz, br, bphi);

    // Print: position, total field, and the first/last basis-mode contributions
    // (a thin probe — the Python driver checks element-wise against mfs).
    std::cout << x << " " << y << " " << z << " "
              << B.x() << " " << B.y() << " " << B.z() << " "
              << bz.front() << " " << br.front() << " " << bphi.front() << " "
              << bz.back()  << " " << br.back()  << " " << bphi.back()  << "\n";
  }

  // Now print the per-mode basis values at one canonical point (so the
  // Python driver can do a full N-mode comparison).
  const double xC = 50.0, yC = 30.0, zC = 100.0;
  const GlobalPoint gpC(xC, yC, zC);
  ev.evaluateBasisAt(gpC, bz, br, bphi);
  std::cout << "# canonical_point " << xC << " " << yC << " " << zC << "\n";
  std::cout << "# canonical_columns: mode_idx l m cs Bz_basis Br_basis Bphi_basis\n";
  for (unsigned int i = 0; i < ev.nModes(); ++i) {
    const auto& p = ev.params()[i];
    std::cout << "@ " << i << " " << p.l << " " << p.m << " " << p.cs << " "
              << bz[i] << " " << br[i] << " " << bphi[i] << "\n";
  }
  return 0;
}
