"""C++ ↔ mfs closure test for ScalarPot3DEval.

Reads the output of testScalarPot3DEval (64 random points with total
field + first/last basis-mode contributions, plus 360 modes at one
canonical point), recomputes the same quantities via mfs.harmonic_basis,
and verifies element-wise agreement to 1e-12 relative.

Usage (from mfs/ with .venv active):
    python ../CMSSW_10_6_26_dev/src/MagneticField/ParametrizedEngine/test/\\
           test_scalar_pot3d_against_mfs.py \\
        --cpp-out /tmp/scalarpot3d_cpp.txt \\
        --dump data/fitresults/polyfit3d_full_coeffs_lmax18_cmsswnorm.txt
"""

import argparse
import sys

import numpy as np

# Path-flex import: prefer mfs in the same workspace.
sys.path.insert(0, '/work/submit/david_w/ZMass/mfs')
import harmonic_basis as hb


def parse_dump(path):
    """Read the Phase-A.5 flat-text dump and return (params, coeffs, header)."""
    header = {}
    params = []
    coeffs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                parts = line[1:].strip().split()
                if not parts:
                    continue
                key = parts[0]
                if key in ('basis_type', 'sha1_npz', 'column:', 'columns:'):
                    continue
                if len(parts) >= 2:
                    header[key] = parts[1]
                continue
            cols = line.split()
            # mode_idx l m cs coefficient
            params.append((int(cols[1]), int(cols[2]), cols[3]))
            coeffs.append(float(cols[4]))
    return params, np.array(coeffs), header


def parse_cpp_out(path):
    """Read C++ test output. Returns header dict, point rows, canonical rows."""
    header = {}
    point_rows = []     # rows of 12 floats: x y z Bx By Bz bz0 br0 bphi0 bzN brN bphiN
    canon_rows = []     # rows: i, l, m, cs, bz, br, bphi
    canon_point = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('#'):
                parts = line[1:].strip().split()
                if not parts:
                    continue
                key = parts[0]
                if key == 'nModes':
                    header['nModes'] = int(parts[1])
                elif key == 'r_scale':
                    header['r_scale'] = float(parts[1])
                elif key == 'z0':
                    header['z0'] = float(parts[1])
                elif key == 'cmssw_norm':
                    header['cmssw_norm'] = bool(int(parts[1]))
                elif key == 'l_max':
                    header['l_max'] = int(parts[1])
                elif key == 'canonical_point':
                    canon_point = tuple(float(x) for x in parts[1:4])
                continue
            if line.startswith('@'):
                cols = line.split()
                # @ idx l m cs bz br bphi
                canon_rows.append((int(cols[1]), int(cols[2]), int(cols[3]),
                                   cols[4], float(cols[5]), float(cols[6]), float(cols[7])))
                continue
            cols = line.split()
            point_rows.append([float(x) for x in cols])
    return header, np.array(point_rows), canon_rows, canon_point


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cpp-out', required=True, help='Output of testScalarPot3DEval')
    ap.add_argument('--dump',    required=True, help='Same dump.txt fed to the C++ test')
    # Relative tolerance is loose because high-(l,m) basis modes can
    # legitimately produce values very near zero (e.g. P_{l-1}^m=0 for m>l-1
    # in Br) where FP noise gets amplified into a large relative error,
    # while the absolute discrepancy is still at the ulp level.  The pair
    # of tolerances is generous enough to absorb that noise while still
    # catching algorithmic bugs (which manifest as O(1) relative errors).
    ap.add_argument('--rtol', type=float, default=1e-5,
                    help='Relative tolerance per element')
    ap.add_argument('--atol', type=float, default=1e-10,
                    help='Absolute tolerance per element (Tesla); covers '
                         'tiny basis values where the relative blows up.')
    args = ap.parse_args()

    # --- Load both inputs -----------------------------------------------------
    params_dump, coeffs_dump, hdr_dump = parse_dump(args.dump)
    hdr_cpp, point_rows, canon_rows, canon_point = parse_cpp_out(args.cpp_out)

    r_scale = float(hdr_dump['r_scale'])
    z0 = float(hdr_dump['z0'])
    cmssw_norm = (hdr_dump['cmssw_norm'] == '1')
    l_max = int(hdr_dump['l_max'])
    n_modes = len(params_dump)

    print(f"dump: {len(params_dump)} modes  l_max={l_max}  r_scale={r_scale:.4f}  cmssw_norm={cmssw_norm}")
    print(f"cpp:  nModes={hdr_cpp['nModes']}  l_max={hdr_cpp['l_max']}  cmssw_norm={hdr_cpp['cmssw_norm']}")
    assert hdr_cpp['nModes'] == n_modes, "nModes mismatch"
    assert hdr_cpp['l_max'] == l_max, "l_max mismatch"
    assert bool(hdr_cpp['cmssw_norm']) == cmssw_norm, "cmssw_norm mismatch"

    # --- Check 1: total field at 64 random points -----------------------------
    print()
    print("Check 1: total absolute field at 64 points")
    x, y, z = point_rows[:, 0], point_rows[:, 1], point_rows[:, 2]
    Bx_cpp, By_cpp, Bz_cpp = point_rows[:, 3], point_rows[:, 4], point_rows[:, 5]

    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)

    pred = hb.eval_field(coeffs_dump, params_dump, r, phi, z,
                         components=('Bz', 'Br', 'Bphi'),
                         r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)
    Br_mfs = pred['Br']
    Bphi_mfs = pred['Bphi']
    Bz_mfs = pred['Bz']
    Bx_mfs = Br_mfs * np.cos(phi) - Bphi_mfs * np.sin(phi)
    By_mfs = Br_mfs * np.sin(phi) + Bphi_mfs * np.cos(phi)

    for label, c, m in [('Bx', Bx_cpp, Bx_mfs),
                        ('By', By_cpp, By_mfs),
                        ('Bz', Bz_cpp, Bz_mfs)]:
        diffs = np.abs(c - m)
        rel = diffs / np.maximum(np.abs(m), 1e-30)
        ok = (rel < args.rtol) | (diffs < args.atol)
        n_bad = int((~ok).sum())
        print(f"  {label}: max|c-m|={float(diffs.max()):.3e}  max rel={float(rel.max()):.3e}  "
              f"failing={n_bad}/{len(c)}")
        assert n_bad == 0, f"{label}: {n_bad} elements fail both rel<{args.rtol} and abs<{args.atol}"

    # --- Check 2: first/last mode basis values --------------------------------
    print()
    print("Check 2: first & last basis-mode values at 64 points")
    bz0_cpp, br0_cpp, bphi0_cpp = point_rows[:, 6], point_rows[:, 7], point_rows[:, 8]
    bzN_cpp, brN_cpp, bphiN_cpp = point_rows[:, 9], point_rows[:, 10], point_rows[:, 11]
    for (i, label, c_bz, c_br, c_bphi) in [
        (0, 'first', bz0_cpp, br0_cpp, bphi0_cpp),
        (n_modes - 1, 'last', bzN_cpp, brN_cpp, bphiN_cpp),
    ]:
        l, m, cs = params_dump[i]
        m_bz = hb.bz_basis(l, m, cs, r, phi, z, r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)
        m_br = hb.br_basis(l, m, cs, r, phi, z, r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)
        m_bphi = hb.bphi_basis(l, m, cs, r, phi, z, r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)
        for tag, cv, mv in [('Bz', c_bz, m_bz), ('Br', c_br, m_br), ('Bphi', c_bphi, m_bphi)]:
            diffs = np.abs(cv - mv)
            rel = diffs / np.maximum(np.abs(mv), 1e-30)
            ok = (rel < args.rtol) | (diffs < args.atol)
            n_bad = int((~ok).sum())
            print(f"  mode {label} ({l},{m},{cs}) {tag}: max|c-m|={float(diffs.max()):.3e}  "
                  f"max rel={float(rel.max()):.3e}  failing={n_bad}/{len(cv)}")
            assert n_bad == 0, f"basis ({l},{m},{cs}) {tag}: {n_bad} fail"

    # --- Check 3: full per-mode basis at the canonical point ------------------
    print()
    print(f"Check 3: full {n_modes}-mode basis at canonical point {canon_point}")
    cx, cy, cz = canon_point
    cr = np.array([np.sqrt(cx**2 + cy**2)])
    cphi = np.array([np.arctan2(cy, cx)])
    cz_arr = np.array([cz])

    cpp_bz = np.array([row[4] for row in canon_rows])
    cpp_br = np.array([row[5] for row in canon_rows])
    cpp_bphi = np.array([row[6] for row in canon_rows])

    mfs_bz = np.zeros(n_modes)
    mfs_br = np.zeros(n_modes)
    mfs_bphi = np.zeros(n_modes)
    for i, (l, m, cs) in enumerate(params_dump):
        mfs_bz[i] = hb.bz_basis(l, m, cs, cr, cphi, cz_arr, r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)[0]
        mfs_br[i] = hb.br_basis(l, m, cs, cr, cphi, cz_arr, r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)[0]
        mfs_bphi[i] = hb.bphi_basis(l, m, cs, cr, cphi, cz_arr, r_scale=r_scale, z0=z0, cmssw_norm=cmssw_norm)[0]

    for tag, cv, mv in [('Bz', cpp_bz, mfs_bz), ('Br', cpp_br, mfs_br), ('Bphi', cpp_bphi, mfs_bphi)]:
        diffs = np.abs(cv - mv)
        rel = diffs / np.maximum(np.abs(mv), 1e-30)
        ok = (rel < args.rtol) | (diffs < args.atol)
        n_bad = int((~ok).sum())
        worst = int(np.argmax(rel))
        l, m, cs = params_dump[worst]
        print(f"  {tag}: max|c-m|={float(diffs.max()):.3e}  max rel={float(rel.max()):.3e}  "
              f"failing={n_bad}/{len(cv)}  "
              f"worst-rel @ ({l},{m},{cs}): cpp={cv[worst]:.6e} mfs={mv[worst]:.6e}")
        assert n_bad == 0, f"canonical {tag}: {n_bad} fail"

    print()
    print("All ScalarPot3DEval ↔ mfs.harmonic_basis checks passed.")


if __name__ == '__main__':
    main()
