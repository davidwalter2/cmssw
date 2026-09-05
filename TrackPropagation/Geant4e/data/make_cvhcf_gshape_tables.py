#!/usr/bin/env python3
"""Dump the offline Moliere shape tables that `cvhcf` (CvhCfExponents) reads.

WHY A DATA FILE AND NOT A RUNTIME BUILD.  `cf_ms_exact._build_elec_tables()`
is 141 rows x 1600 tau x (2048 + 512) quadrature nodes = 5.8e8 evaluations of
J0(x) - 1.  Rebuilding it in C++ costs ~20 s of startup AND makes the answer
depend on which J0 implementation the compiler ships -- libstdc++'s
`std::cyl_bessel_j` is not Cephes, which is what scipy uses, so the two tables
would differ by more than the 1e-6 the in-maker exponent is required to hold
against the offline reference.  Dumping the reference's OWN table removes both
problems: the C++ evaluates the identical numbers with the identical linear
interpolation, and the only remaining difference is the arithmetic of the sum.

ONLY THE TWO KERNELS THE PRODUCTION PATH USES are written:
  `dipole` -- the NUCLEAR term, because `cf_track_resolution.MS_FINE_G = 1.0`
              routes the nuclear piece through `gshape_elec(..., "dipole")`
              rather than through `gshape`;
  `hard`   -- the ATOMIC-ELECTRON term at its own ceiling, because
              `MS_ELEC_TMAX = 1.0` and `MS_ELEC_EDGE = 1.0`.
`kine0`/`kine1` are diagnostics (MS_ELEC_EDGE 2/3) and are NOT written; if a
study ever needs them, add them here and bump the format version.

The delta-ray recoil table `cf_delta_ray._PHI` (4001 nodes, 32 kB) rides in the
same file for the same reason: it too is a cumulative quadrature of J0(x) - 1.

usage:
    source /work/submit/david_w/ZMass/mfs/.venv/bin/activate
    python3 make_cvhcf_gshape_tables.py [-o cvhcf_gshape_elec_v1.bin]
"""
import argparse
import hashlib
import os
import struct
import sys

import numpy as np

RES = "/work/submit/david_w/ZMass/calibration_studies/resolution"
MAGIC = b"CVHCFGSH"
VERSION = 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "cvhcf_gshape_elec_v1.bin"))
    a = p.parse_args()

    sys.path.insert(0, RES)
    import cf_ms_exact as M                                  # noqa: E402
    import cf_delta_ray as D                                 # noqa: E402

    gtau = np.asarray(M._GTAU, dtype=np.float64)
    ely = np.asarray(M._ELEC_Y, dtype=np.float64)
    hard = np.ascontiguousarray(M._GE["hard"], dtype=np.float64)
    dip = np.ascontiguousarray(M._GE["dipole"], dtype=np.float64)
    assert hard.shape == (len(ely), len(gtau)) and dip.shape == hard.shape

    plx = np.asarray(D._LX, dtype=np.float64)
    pphi = np.asarray(D._PHI, dtype=np.float64)
    assert plx.shape == pphi.shape

    with open(a.out, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<iiii", VERSION, len(gtau), len(ely), len(plx)))
        f.write(gtau.tobytes())
        f.write(ely.tobytes())
        f.write(hard.tobytes())
        f.write(dip.tobytes())
        f.write(plx.tobytes())
        f.write(pphi.tobytes())
    n = os.path.getsize(a.out)
    md5 = hashlib.md5(open(a.out, "rb").read()).hexdigest()
    print(f"wrote {a.out}  {n} bytes  md5 {md5}")
    print(f"  gtau  {len(gtau)}  [{gtau[0]:g}, {gtau[-1]:g}]")
    print(f"  elecY {len(ely)}  [{ely[0]:g}, {ely[-1]:g}]")
    print(f"  philx {len(plx)}  [{plx[0]:g}, {plx[-1]:g}]  "
          f"Phi[0]={pphi[0]:.17g} Phi[-1]={pphi[-1]:.17g}")
    print(f"  J0M1_GUARD={M.J0M1_GUARD} G4_FF_SQUARED={M.G4_FF_SQUARED} "
          f"MS_CHI0_G4={M.MS_CHI0_G4} MS_FF_G4={M.MS_FF_G4}")


if __name__ == "__main__":
    main()
