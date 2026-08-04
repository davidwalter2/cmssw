#!/usr/bin/env python3
"""Material grouping v2, built on three axes (2026-07-31 design):

1. LOCATION -- the correction must act at the right point along the
   trajectory (V0 daughters start displaced, cosmics run top-to-bottom,
   Z muons come from the IP). Implemented by keeping the validated
   per-layer/per-ring active skeleton of materialGroups50 and keeping all
   passive groups subsystem-local (with r/z windows where a subsystem
   spans distinct regions).

2. MATERIAL -- one parameter must describe energy loss AND multiple
   scattering coherently, which requires composition homogeneity
   (Z/A and Z(Z+1)/A fixed within a group). A Geant4 logical volume has
   exactly one material, so volume-name rules give EXACT material purity.
   Each passive location group is split into physics classes:
     composite (CF/polymer/Nomex/epoxy, X0_mass 30-46 g/cm2)
     metal     (Al/Cu/steel-dominated, X0_mass < 30)
     coolant   (C6F14/water/silicone, name-matched)
     air       (known gas)

3. UNCERTAINTY -- float what is uncertain, pin what is known:
     air            prior 0.005  (rho*d from T/P; effectively fixed)
     active silicon prior 0.02   (sensor thickness known)
     beampipe (Be)  prior 0.02
     metal          prior 0.08   (machined/drawn parts)
     coolant        prior 0.10
     composite      prior 0.20   (layup/adhesive budgets -- the dominant
                                  uncertainty class, 7.7% of dE in CF alone)
     misc/other     prior 0.10

Usage: makeMaterialGroupsV2.py <material_audit.txt> <outfile>
"""
import re
import sys
from collections import defaultdict

PRIORS = {"air": 0.005, "active": 0.02, "beampipe": 0.02, "metal": 0.08,
          "coolant": 0.10, "composite": 0.20, "misc": 0.10}

# location skeleton kept from materialGroups50 (validated): actives per
# layer/ring + beampipe. (name, regex, rmin, rmax, zmin, zmax, prior)
SKELETON = [
    # Beryllium WALL only: the BEAM/BeamVacuum volumes inside are air and
    # vacuum, and a wide regex here swallows them before the air rule can
    # claim them (found 2026-08-02: the group then inherits the impure-group
    # artifact, -0.55, where the pure Be wall closes at -0.19+-0.12).
    ("beampipe",       r"^BeamTube", None, 3.5, None, None, PRIORS["beampipe"]),
    ("bpix_active_L1", r"^PixelBarrel(Active|Sensor|ROChip|TBM)", None, 5.9, None, 30., PRIORS["active"]),
    ("bpix_active_L2", r"^PixelBarrel(Active|Sensor|ROChip|TBM)", 5.9, 8.7, None, 30., PRIORS["active"]),
    ("bpix_active_L3", r"^PixelBarrel(Active|Sensor|ROChip|TBM)", 8.7, None, None, 30., PRIORS["active"]),
    ("fpix_active_D1", r"^PixelForwardActive", None, None, None, 38., PRIORS["active"]),
    ("fpix_active_D2", r"^PixelForwardActive", None, None, 38., None, PRIORS["active"]),
    ("tib_active_L1",  r"^TIB(Active|Hybrid)", None, 30., None, None, PRIORS["active"]),
    ("tib_active_L2",  r"^TIB(Active|Hybrid)", 30., 38., None, None, PRIORS["active"]),
    ("tib_active_L3",  r"^TIB(Active|Hybrid)", 38., 46., None, None, PRIORS["active"]),
    ("tib_active_L4",  r"^TIB(Active|Hybrid)", 46., None, None, None, PRIORS["active"]),
    ("tid_active_R1",  r"^TIDModule0.*Active", None, None, None, None, PRIORS["active"]),
    ("tid_active_R2",  r"^TIDModule1.*Active", None, None, None, None, PRIORS["active"]),
    ("tid_active_R3",  r"^TIDModule2.*Active", None, None, None, None, PRIORS["active"]),
    ("tob_active_L1",  r"^TOBActive", None, 65., None, None, PRIORS["active"]),
    ("tob_active_L2",  r"^TOBActive", 65., 73., None, None, PRIORS["active"]),
    ("tob_active_L3",  r"^TOBActive", 73., 82., None, None, PRIORS["active"]),
    ("tob_active_L4",  r"^TOBActive", 82., 92., None, None, PRIORS["active"]),
    ("tob_active_L5",  r"^TOBActive", 92., 102., None, None, PRIORS["active"]),
    ("tob_active_L6",  r"^TOBActive", 102., None, None, None, PRIORS["active"]),
    ("tec_active_R1",  r"^TECModule0.*Active", None, None, None, None, PRIORS["active"]),
    ("tec_active_R2",  r"^TECModule1.*Active", None, None, None, None, PRIORS["active"]),
    ("tec_active_R3",  r"^TECModule2.*Active", None, None, None, None, PRIORS["active"]),
    ("tec_active_R4",  r"^TECModule3.*Active", None, None, None, None, PRIORS["active"]),
    ("tec_active_R5",  r"^TECModule4.*Active", None, None, None, None, PRIORS["active"]),
    ("tec_active_R6",  r"^TECModule5.*Active", None, None, None, None, PRIORS["active"]),
    ("tec_active_R7",  r"^TECModule6.*Active", None, None, None, None, PRIORS["active"]),
]
SKEL_RE = [re.compile(r[1]) for r in SKELETON]

COOLANT_RE = re.compile(r"[Cc]oolant|C6F14|Water|silicone_pipes|Silicone_Gel")
AIR_RE = re.compile(r"^(T_)?(E_)?Air$")

# passive location regions: subsystem prefix from the LV name
REGIONS = [
    ("bpix",   re.compile(r"^(PixelBarrel|PIXPatchpanel)")),
    ("fpix",   re.compile(r"^(PixelForward|PixelFwd|TPGBlade|CarbonFiberSki)")),
    ("tibtid", re.compile(r"^TIBTID")),
    ("tib",    re.compile(r"^TIB")),
    ("tid",    re.compile(r"^TID")),
    ("tob",    re.compile(r"^TOB")),
    ("tec",    re.compile(r"^(TEC|BHDisk|BHCovers)")),
    ("tracker", re.compile(r"")),   # catch-all: Tracker, support tube, thermal screen, SF*, EI*
]


def matclass(material, x0mass):
    if AIR_RE.search(material):
        return "air"
    if COOLANT_RE.search(material):
        return "coolant"
    if x0mass < 30.:
        return "metal"
    if x0mass < 46.:
        return "composite"
    return "misc"


def main(auditfile, outfile):
    # volume -> (material, dE, mass, x0) aggregated over bins
    vol = defaultdict(lambda: [None, 0., 0., 0.])
    for line in open(auditfile):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        name, mat = f[0], f[9]
        x0, mass, dE = float(f[4]), float(f[5]), float(f[6])
        v = vol[name]
        v[0] = mat
        v[1] += dE
        v[2] += mass
        v[3] += x0

    # drop what the skeleton already covers
    passive = {}
    for name, (mat, dE, mass, x0) in vol.items():
        if any(rx.search(name) for rx in SKEL_RE):
            continue
        x0mass = mass / x0 if x0 > 0 else 40.
        passive[name] = (mat, matclass(mat, x0mass), dE, x0)

    # region x class aggregation
    groups = defaultdict(lambda: [[], 0., 0.])   # names, dE, x0
    for name, (mat, cls, dE, x0) in passive.items():
        for rname, rrx in REGIONS:
            if rrx.search(name):
                break
        # air and coolant are global classes (one group each) -- their
        # physics and uncertainty do not depend on which subsystem they
        # sit in, and per-region air groups would be unconstrained
        key = cls if cls in ("air", "coolant") else f"{rname}_{cls}"
        g = groups[key]
        g[0].append(name)
        g[1] += dE
        g[2] += x0

    tot_dE = sum(v[1] for v in vol.values())
    # merge groups below threshold into the region misc
    MIN_FRAC = 0.002
    final = {}
    for key, (names, dE, x0) in sorted(groups.items(), key=lambda kv: -kv[1][1]):
        if key in ("air", "coolant") or dE / tot_dE >= MIN_FRAC:
            final[key] = (names, dE, x0)
        else:
            region = key.split("_")[0]
            misc = f"{region}_misc"
            if misc not in final:
                final[misc] = ([], 0., 0.)
            final[misc][0].extend(names)
            final[misc] = (final[misc][0], final[misc][1] + dE, final[misc][2] + x0)

    with open(outfile, "w") as out:
        out.write("# materialGroupsV2 -- location x material-class x uncertainty-class\n")
        out.write(f"# generated by makeMaterialGroupsV2.py from {auditfile}\n")
        out.write("# RULE\tgroupId\tgroupName\tnameRegex\trmin\trmax\tzmin\tzmax\tzside\tk_init\tprior_sigma\n")
        out.write("# first matching rule wins; unmatched volumes -> group 0 (other)\n")
        gid = 1

        def wr(name, regex, rmin, rmax, zmin, zmax, prior):
            nonlocal gid
            fmt = lambda v: "-" if v is None else f"{v:g}"
            out.write(f"RULE\t{gid}\t{name}\t{regex}\t{fmt(rmin)}\t{fmt(rmax)}"
                      f"\t{fmt(zmin)}\t{fmt(zmax)}\t0\t0.0\t{prior}\n")
            gid += 1

        for r in SKELETON:
            wr(r[0], r[1], r[2], r[3], r[4], r[5], r[6])
        # class groups: exact volume-name alternation = exact material purity
        for key, (names, dE, x0) in sorted(final.items(), key=lambda kv: -kv[1][1]):
            cls = key.split("_")[-1] if key not in ("air", "coolant") else key
            prior = PRIORS.get(cls, PRIORS["misc"])
            regex = "^(" + "|".join(sorted(re.escape(n) for n in set(names))) + ")$"
            wr(key, regex, None, None, None, None, prior)
            print(f"{key:22s} dE {100*dE/tot_dE:5.2f}%  x/X0-share {x0:8.1f}  "
                  f"{len(set(names)):3d} volumes  prior {prior}")
    print(f"\nwrote {outfile}: {gid-1} rules "
          f"({len(SKELETON)} skeleton + {gid-1-len(SKELETON)} class groups)")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
