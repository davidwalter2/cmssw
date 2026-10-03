#!/usr/bin/env python3
"""Dump the offline Moliere shape tables that `cvhcf` (CvhCfExponents) reads.

WHY A DATA FILE AND NOT A RUNTIME BUILD.  `cf_ms_exact._build_elec_tables()`
is 171 rows x 1600 tau x (2048 + 512) quadrature nodes = 7e8 evaluations of
J0(x) - 1.  Rebuilding it in C++ costs ~25 s of startup AND makes the answer
depend on which J0 implementation the compiler ships.  Dumping the
reference's OWN table removes both problems: the C++ evaluates the identical
numbers with the identical linear interpolation, and the only remaining
difference is the arithmetic of the sum.

ONLY THE TWO KERNELS THE PRODUCTION PATH USES are written:
  `dipole` -- the NUCLEAR term, because `cf_track_resolution.MS_FINE_G = 1.0`
              routes the nuclear piece through `gshape_elec(..., "dipole")`
              rather than through `gshape`;
  `hard`   -- the ATOMIC-ELECTRON term at its own ceiling, because
              `MS_ELEC_TMAX = 1.0` and `MS_ELEC_EDGE = 1.0`.
`kine0`/`kine1` are diagnostics (MS_ELEC_EDGE 2/3) and are NOT written.

The rows run over `cf_ms_exact._ELEC_Y` = logspace(-1, 16, 171) in
Y = ymax^2: from y = 10^-1/2, where the electron ceiling sits when the
knock-on channel carries the collisions above the e- production threshold
(cf_knockon.KNOCKON_JOINT), up to the kinematic ceiling of any momentum.

BINARY LAYOUT (little-endian): char[8] "CVHCFGSH", i4 version = 2,
i4 nTau, i4 nY, f8[nTau] gtau, f8[nY] elecY, f8[nY][nTau] hard,
f8[nY][nTau] dipole.  The file name stem is the table id `modelTag` records.

usage:
    source /work/submit/david_w/ZMass/calibration_studies/setup_env.sh
    python3 make_cvhcf_gshape_tables.py [-o cvhcf_gshape_elec_v2.bin]
"""
import argparse
import hashlib
import os
import struct
import sys

import numpy as np

RES = "/work/submit/david_w/ZMass/calibration_studies/resolution"
MAGIC = b"CVHCFGSH"
VERSION = 2


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "cvhcf_gshape_elec_v2.bin"))
    a = p.parse_args()

    sys.path.insert(0, RES)
    import cf_ms_exact as M                                  # noqa: E402

    gtau = np.asarray(M._GTAU, dtype=np.float64)
    ely = np.asarray(M._ELEC_Y, dtype=np.float64)
    hard = np.ascontiguousarray(M._GE["hard"], dtype=np.float64)
    dip = np.ascontiguousarray(M._GE["dipole"], dtype=np.float64)
    assert hard.shape == (len(ely), len(gtau)) and dip.shape == hard.shape

    with open(a.out, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<iii", VERSION, len(gtau), len(ely)))
        f.write(gtau.tobytes())
        f.write(ely.tobytes())
        f.write(hard.tobytes())
        f.write(dip.tobytes())
    n = os.path.getsize(a.out)
    md5 = hashlib.md5(open(a.out, "rb").read()).hexdigest()
    print(f"wrote {a.out}  {n} bytes  md5 {md5}")
    print(f"  gtau  {len(gtau)}  [{gtau[0]:g}, {gtau[-1]:g}]")
    print(f"  elecY {len(ely)}  [{ely[0]:g}, {ely[-1]:g}]")
    print(f"  J0M1_GUARD={M.J0M1_GUARD} G4_FF_SQUARED={M.G4_FF_SQUARED} "
          f"MS_CHI0_G4={M.MS_CHI0_G4} MS_FF_G4={M.MS_FF_G4}")


if __name__ == "__main__":
    main()
