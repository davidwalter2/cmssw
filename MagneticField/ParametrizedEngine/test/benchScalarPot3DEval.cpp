// Microbenchmark: time evaluateAbsoluteAt over 100k random points.
// Build via the test BuildFile.xml; run inside el7 container.
//
//   benchScalarPot3DEval <dump.txt>

#include <chrono>
#include <cstdio>
#include <random>
#include <vector>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/ParametrizedEngine/interface/ScalarPot3DEval.h"

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::fprintf(stderr, "usage: %s <dump.txt>\n", argv[0]);
    return 1;
  }
  magfieldparam::ScalarPot3DEval ev(argv[1]);
  const auto& c = ev.initCoeffs();

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

  // Warmup: 1k iterations to let TLBs/caches settle.
  double sink = 0.0;
  for (unsigned int k = 0; k < 1000; ++k) {
    auto B = ev.evaluateAbsoluteAt(pts[k % N_PTS], c);
    sink += B.x() + B.y() + B.z();
  }

  // Hot loop.
  auto t0 = std::chrono::steady_clock::now();
  for (unsigned int k = 0; k < N_PTS; ++k) {
    auto B = ev.evaluateAbsoluteAt(pts[k], c);
    sink += B.x() + B.y() + B.z();
  }
  auto t1 = std::chrono::steady_clock::now();
  const double secs = std::chrono::duration<double>(t1 - t0).count();
  const double us_per_call = secs * 1e6 / static_cast<double>(N_PTS);

  std::printf("evaluateAbsoluteAt: %u calls in %.3f s -> %.2f us/call  (sink=%.6e)\n",
              N_PTS, secs, us_per_call, sink);
  return 0;
}
