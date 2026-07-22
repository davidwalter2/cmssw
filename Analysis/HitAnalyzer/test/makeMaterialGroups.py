#!/usr/bin/env python3
"""Phase 0 of the global material model: build the two grouping tiers
(materialGroups50 / materialGroups100) from the MaterialAuditAnalyzer
per-volume tally.

Groups are defined as ordered classification rules
    (name regex, r window [cm], |z| window [cm], zside)  ->  groupId
because logical-volume names alone cannot resolve layers (the same LV is
placed at several layers, e.g. PixelBarrelActive* at all three BPIX
layers) -- the runtime classifier in the propagator applies the same
rules per step using the step midpoint.

Tier 50: z-symmetric, phi-symmetric (zside=0 everywhere).
Tier 100: the same taxonomy with z+/z- split for everything except the
beam pipe (endcap elements, services, and barrel groups alike).

Output format (one file per tier):
  RULE  <groupId>  <groupName>  <nameRegex>  <rmin> <rmax>  <zmin> <zmax>  <zside>  <k_init>  <prior_sigma>
First matching rule wins; group 0 = catch-all "other".

Usage: makeMaterialGroups.py <material_audit.txt> <outdir>
"""
import re
import sys
from collections import defaultdict

# taxonomy: (name, regex, rmin, rmax, zmin, zmax, prior)
# r/z in cm on the step midpoint; None = unbounded. Ordered: first match wins.
BASE_RULES = [
    ("beampipe",        r"^(BeamTube|BeamVacuum|BEAM)", None, 3.5, None, None, 0.01),
    # --- pixel barrel ---
    ("bpix_active_L1",  r"^PixelBarrel(Active|Sensor|ROChip|TBM)", None, 5.9, None, 30., 0.02),
    ("bpix_active_L2",  r"^PixelBarrel(Active|Sensor|ROChip|TBM)", 5.9, 8.7, None, 30., 0.02),
    ("bpix_active_L3",  r"^PixelBarrel(Active|Sensor|ROChip|TBM)", 8.7, None, None, 30., 0.02),
    ("bpix_services",   r"^(PixelBarrel(Conn|Supply|Flange)|PIXPatchpanel)", None, None, None, None, 0.10),
    ("bpix_support",    r"^PixelBarrel", None, None, None, 30., 0.05),
    ("bpix_support",    r"^PixelBarrel", None, None, 30., None, 0.10),
    # --- pixel forward ---
    ("fpix_active_D1",  r"^PixelForwardActive", None, None, None, 38., 0.02),
    ("fpix_active_D2",  r"^PixelForwardActive", None, None, 38., None, 0.02),
    ("fpix_support",    r"^(PixelForward|PixelFwd)", None, None, None, None, 0.05),
    # --- TIB (layer split by r: L1 ~25, L2 ~34, L3 ~42, L4 ~50) ---
    ("tib_active_L1",   r"^TIB(Active|Hybrid)", None, 30., None, None, 0.02),
    ("tib_active_L2",   r"^TIB(Active|Hybrid)", 30., 38., None, None, 0.02),
    ("tib_active_L3",   r"^TIB(Active|Hybrid)", 38., 46., None, None, 0.02),
    ("tib_active_L4",   r"^TIB(Active|Hybrid)", 46., None, None, None, 0.02),
    ("tibtid_services", r"^TIBTID", None, None, None, None, 0.10),
    ("tib_support",     r"^TIB", None, None, None, None, 0.05),
    # --- TID (ring encoded in the name) ---
    ("tid_active_R1",   r"^TIDModule0.*Active", None, None, None, None, 0.02),
    ("tid_active_R2",   r"^TIDModule1.*Active", None, None, None, None, 0.02),
    ("tid_active_R3",   r"^TIDModule2.*Active", None, None, None, None, 0.02),
    ("tid_support",     r"^TID", None, None, None, None, 0.05),
    # --- TOB (name digits are module TYPES spanning layer pairs -- the
    #     same LV serves two layers; split layers by r instead:
    #     L1..L6 at r = 61, 69, 78, 87, 97, 108 cm) ---
    ("tob_active_L1",   r"^TOBActive", None, 65., None, None, 0.02),
    ("tob_active_L2",   r"^TOBActive", 65., 73., None, None, 0.02),
    ("tob_active_L3",   r"^TOBActive", 73., 82., None, None, 0.02),
    ("tob_active_L4",   r"^TOBActive", 82., 92., None, None, 0.02),
    ("tob_active_L5",   r"^TOBActive", 92., 102., None, None, 0.02),
    ("tob_active_L6",   r"^TOBActive", 102., None, None, None, 0.02),
    ("tob_services",    r"^TOB(Cable|AxServices|RadServices)", None, None, None, None, 0.10),
    ("tob_support",     r"^TOB", None, None, None, None, 0.05),
    # --- TEC (ring encoded in the name for actives) ---
    ("tec_active_R1",   r"^TECModule0.*Active", None, None, None, None, 0.02),
    ("tec_active_R2",   r"^TECModule1.*Active", None, None, None, None, 0.02),
    ("tec_active_R3",   r"^TECModule2.*Active", None, None, None, None, 0.02),
    ("tec_active_R4",   r"^TECModule3.*Active", None, None, None, None, 0.02),
    ("tec_active_R5",   r"^TECModule4.*Active", None, None, None, None, 0.02),
    ("tec_active_R6",   r"^TECModule5.*Active", None, None, None, None, 0.02),
    ("tec_active_R7",   r"^TECModule6.*Active", None, None, None, None, 0.02),
    ("tec_services",    r"^TEC(Cool|ServChannel|PpBox|PpConnector)", None, None, None, None, 0.10),
    ("tec_structure",   r"^TEC", None, None, None, None, 0.05),
    # --- common structures ---
    ("thermal_screen",  r"^TrackerThe", None, None, None, None, 0.05),
    ("support_tube",    r"^TrackerSup", None, None, None, None, 0.05),
    ("pp1_cables",      r"^Tracker_PP1", None, None, None, None, 0.10),
    ("pixel_patch",     r"^(MCHead|PixelServ)", None, None, None, None, 0.10),
]
OTHER = ("other", None, None, None, None, None, 0.20)


def classify(name, r, z):
    for i, (gname, rex, rmin, rmax, zmin, zmax, prior) in enumerate(BASE_RULES):
        if not re.match(rex, name):
            continue
        if rmin is not None and r < rmin:
            continue
        if rmax is not None and r >= rmax:
            continue
        if zmin is not None and abs(z) < zmin:
            continue
        if zmax is not None and abs(z) >= zmax:
            continue
        return i
    return -1  # other


def main():
    audit, outdir = sys.argv[1], sys.argv[2]
    rows = []
    for line in open(audit):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        rows.append((f[0], int(f[1]), float(f[3]), float(f[4]), float(f[6]),
                     float(f[7]), float(f[8])))  # name side path xX0 dE meanR meanZ

    # tier 50: zside-merged groups; tier 100: everything except the beam
    # pipe split by zside
    for tier, zsplit in (("50", False), ("100", True)):
        gsum = defaultdict(lambda: [0., 0.])  # (ruleIdx, zside) -> [xX0, dE]
        for name, side, path, xX0, dE, r, z in rows:
            idx = classify(name, r, z)
            key = (idx, side if (zsplit and idx >= 0) else 0)
            gsum[key][0] += xX0
            gsum[key][1] += dE

        # assign contiguous group ids: 0 = other, then rule order (z- before z+)
        keys = sorted(k for k in gsum if k[0] >= 0)
        gid = {}
        names = {}
        nextid = 1
        for k in keys:
            gid[k] = nextid
            rule = BASE_RULES[k[0]]
            suffix = {0: "", -1: "_zm", 1: "_zp"}[k[1]]
            names[k] = rule[0] + suffix
            nextid += 1

        out = f"{outdir}/materialGroups{tier}.txt"
        with open(out, "w") as fo:
            fo.write(f"# materialGroups{tier} v1 -- generated by makeMaterialGroups.py from {audit}\n")
            fo.write("# RULE  groupId  groupName  nameRegex  rmin rmax  zmin zmax  zside  k_init  prior_sigma\n")
            fo.write("# first matching rule wins; unmatched volumes -> group 0 (other)\n")
            for k in keys:
                i, side = k
                gname, rex, rmin, rmax, zmin, zmax, prior = BASE_RULES[i]
                fmt = lambda v: "-" if v is None else f"{v:g}"
                fo.write(f"RULE\t{gid[k]}\t{names[k]}\t{rex}\t{fmt(rmin)}\t{fmt(rmax)}\t"
                         f"{fmt(zmin)}\t{fmt(zmax)}\t{side}\t0.0\t{prior:g}\n")

        # summary merged by group name (several rules can share a name)
        bysum = defaultdict(lambda: [0., 0.])
        for k, v in gsum.items():
            bysum[names.get(k, "other")][0] += v[0]
            bysum[names.get(k, "other")][1] += v[1]
        total = sum(v[0] for v in bysum.values())
        totdE = sum(v[1] for v in bysum.values())
        print(f"=== tier {tier}: {nextid - 1} groups + other -> {out}")
        print(f"    {'group':24s} {'x/X0 %':>8s} {'dE %':>8s}")
        for nm in sorted(bysum, key=lambda n: -bysum[n][0]):
            print(f"    {nm:24s} {100*bysum[nm][0]/total:8.2f} {100*bysum[nm][1]/totdE:8.2f}")


if __name__ == "__main__":
    main()
