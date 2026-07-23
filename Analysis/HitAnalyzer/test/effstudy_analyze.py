#!/usr/bin/env python3
"""
A/B efficiency analysis for the eb96caef glued-module fix.

Reads the two single-track CVH-refit outputs produced by
runCvhSingleTrackEffStudy.py on the SAME events:
  effstudy_fix_0.root   (gluedTiltThr=0.05 -> repair ON)
  effstudy_bug_0.root   (gluedTiltThr=1e9  -> repair OFF, buggy)

Since a failed fit writes no tree row and both arms see identical input,
the per-bin denominator is identical and cancels, so

    SF_bug(bin) = eps_buggy / eps_fixed = N_pass_bug(bin) / N_pass_fix(bin)

is applied as a per-muon weight to MC. This also matches tracks that pass
BOTH arms to read off the momentum bias the fix induces on survivors
(the check an efficiency SF alone cannot capture).

Usage:
  python effstudy_analyze.py --fix effstudy_fix_0.root --bug effstudy_bug_0.root \
      [--ptmin 25 --ptmax 65]
"""
import argparse
import numpy as np
import uproot

BR = ["run", "lumi", "event", "trackPt", "trackEta", "trackPhi",
      "trackCharge", "edmvalref", "niter"]


def load(paths):
    # paths: one or more files/globs; concatenate the 'tree' across all.
    import glob
    files = []
    for p in paths:
        files += sorted(glob.glob(p)) if any(c in p for c in "*?[") else [p]
    if not files:
        raise SystemExit(f"no files matched: {paths}")
    d = uproot.concatenate([f"{f}:tree" for f in files],
                           BR + ["refParms"], library="np")
    d["qop"] = np.array([rp[0] for rp in d["refParms"]], dtype=float)
    print(f"  loaded {len(files)} file(s), {len(d['run'])} tracks")
    return d


def key(d, i):
    # stable per-track id: event coordinates + input (run-independent) kinematics
    return (int(d["run"][i]), int(d["lumi"][i]), int(d["event"][i]),
            round(float(d["trackEta"][i]), 3), round(float(d["trackPhi"][i]), 3))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fix", required=True, nargs="+", help="fix file(s)/glob(s)")
    ap.add_argument("--bug", required=True, nargs="+", help="bug file(s)/glob(s)")
    ap.add_argument("--ptmin", type=float, default=0.0)
    ap.add_argument("--ptmax", type=float, default=1e9)
    a = ap.parse_args()

    fix, bug = load(a.fix), load(a.bug)

    def sel(d):
        m = (d["trackPt"] >= a.ptmin) & (d["trackPt"] < a.ptmax)
        return {k: (v[m] if hasattr(v, "__len__") and len(v) == len(m) else v)
                for k, v in d.items()}
    fix, bug = sel(fix), sel(bug)

    kf = {key(fix, i): i for i in range(len(fix["run"]))}
    kb = {key(bug, i): i for i in range(len(bug["run"]))}
    only_fix = set(kf) - set(kb)   # pass with fix, fail buggy -> recovered by fix
    only_bug = set(kb) - set(kf)   # pass buggy, fail with fix
    both = set(kf) & set(kb)

    print(f"pT window [{a.ptmin}, {a.ptmax}) GeV")
    print(f"  passing tracks:  fix={len(kf)}  bug={len(kb)}")
    print(f"  pass fix & fail bug (fix recovers): {len(only_fix)}")
    print(f"  pass bug & fail fix              : {len(only_bug)}")
    print(f"  pass both                        : {len(both)}")
    print(f"  global efficiency ratio bug/fix  : {len(kb)/max(len(kf),1):.5f}")

    if only_fix or only_bug:
        print("\n-- flipped tracks (eta, phi, pt, q) : where the bug bites --")
        for tag, s, d in [("FIXonly", only_fix, fix), ("BUGonly", only_bug, bug)]:
            km = {key(d, i): i for i in range(len(d["run"]))}
            for k in list(s)[:40]:
                i = km[k]
                print(f"  {tag}  eta={d['trackEta'][i]:+.3f} phi={d['trackPhi'][i]:+.3f} "
                      f"pt={d['trackPt'][i]:6.2f} q={int(d['trackCharge'][i]):+d}")

    # scale bias on survivors: dq/p (fix - bug) for tracks passing both
    if both:
        dqop = np.array([fix["qop"][kf[k]] - bug["qop"][kb[k]] for k in both])
        eta = np.array([fix["trackEta"][kf[k]] for k in both])
        phi = np.array([fix["trackPhi"][kf[k]] for k in both])
        print("\n-- survivor scale bias  d(q/p) = fix - bug --")
        print(f"  median |d(q/p)| = {np.median(np.abs(dqop)):.3e}  "
              f"max = {np.max(np.abs(dqop)):.3e}")
        big = np.abs(dqop) > 1e-4
        if big.any():
            print(f"  {big.sum()} survivors with |d(q/p)|>1e-4 (biased-but-kept):")
            for e, p, dq in zip(eta[big], phi[big], dqop[big]):
                print(f"     eta={e:+.3f} phi={p:+.3f} d(q/p)={dq:+.3e}")


if __name__ == "__main__":
    main()
