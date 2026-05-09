// Microbenchmark: time PolyFit3DParametrizedMagneticField::inTesla over
// 100k random points. Apples-to-apples comparison with benchScalarPot3DEval.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/ParametrizedEngine/src/PolyFit3DParametrizedMagneticField.h"

int main() {
  // Match PolyFit3DMf cfi default: 3.81143 T (Run 2 nominal).
  edm::ParameterSet pset;
  pset.addParameter<double>("BValue", 3.81143);
  PolyFit3DParametrizedMagneticField field(pset);

  std::mt19937 rng(20260508);
  std::uniform_real_distribution<double> ur(5.0, 290.0);
  std::uniform_real_distribution<double> uphi(-3.14159265358979, 3.14159265358979);
  std::uniform_real_distribution<double> uz(-300.0, 300.0);

  const unsigned int N_PTS = 100000;
  std::vector<GlobalPoint> pts;
  pts.reserve(N_PTS);
  for (unsigned int k = 0; k < N_PTS; ++k) {
    const double r = ur(rng);
    const double phi = uphi(rng);
    const double z = uz(rng);
    pts.emplace_back(r * std::cos(phi), r * std::sin(phi), z);
  }

  double sink = 0.0;
  for (unsigned int k = 0; k < 1000; ++k) {
    auto B = field.inTesla(pts[k % N_PTS]);
    sink += B.x() + B.y() + B.z();
  }

  auto t0 = std::chrono::steady_clock::now();
  for (unsigned int k = 0; k < N_PTS; ++k) {
    auto B = field.inTesla(pts[k]);
    sink += B.x() + B.y() + B.z();
  }
  auto t1 = std::chrono::steady_clock::now();
  const double secs = std::chrono::duration<double>(t1 - t0).count();
  const double us_per_call = secs * 1e6 / static_cast<double>(N_PTS);

  std::printf("PolyFit3D::inTesla: %u calls in %.3f s -> %.2f us/call  (sink=%.6e)\n",
              N_PTS, secs, us_per_call, sink);
  return 0;
}
