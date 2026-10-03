#!/usr/bin/env python3
"""Build a toy geometry from the REAL material the reference trajectory crosses.

WHY THIS EXISTS
---------------
There are two geometries in the clean-propagation study and they differ in too
many ways at once.  The layered toy is 15 shells of 1 mm at rho = 9 g/cm3 with
VACUUM between them; the real tracker is thin silicon separated by long AIR
gaps, with hybrids, cooling pipes, carbon fibre and cable in between.  Several
open items point at that difference rather than at the physics:

  * the real-geometry q/p closure is +0.0676 at the outermost plane and grows
    monotonically with accumulated material, while the toy's is -0.0011;
  * MS_WVI_SPLIT was gauged on the toy and then hard-failed the real geometry,
    where q^2 U/4 reaches 1.7e5 against a validity ceiling of 60 -- a
    thin-target construction meeting a thin-target-free geometry;
  * the locx plateau is the same size and sign on BOTH geometries, so something
    is common to them.

This builds the missing rung: the real RADIAL MATERIAL SEQUENCE, as coaxial
cylinders.  Same materials, same thicknesses, same step structure, same field --
but no stereo modules, no phi gaps, no module edges, no misalignment.  A closure
difference between this and the layered toy is material arrangement; a
difference between this and the real tracker is one of the things it drops.

HOW
---
The reference propagation already knows the answer: `CVH_DEDX_DEBUG=1` makes
Geant4ePropagator dump, for every transport step,

    [dedxdbg] iter= len= dEdx= mat=<G4 name> rho= r= accum=

with lengths and radii in cm.  This script coalesces consecutive steps of the
SAME material into volumes and emits one cylinder per volume.

COALESCING IS NOT OPTIONAL.  Step boundaries are not geometry boundaries: the
propagator's 10 mm ceiling cuts a long air gap into a dozen steps (the real
pt = 3 traversal has runs of 10 and 12 consecutive Air steps).  Volumes are
where the material NAME changes.

MATERIALS ARE REUSED BY NAME, not by effective Z.  The DDD material sections
are already loaded by the job, so `<rMaterial name="trackermaterial:Silicon"/>`
costs nothing, and effective Z is not good enough for the simulation side: the
nuclear-elastic kernel is per-element and a hydrogen kernel is 3x wider than
oxygen at the same rate, which is a 10x over-correction when 8.8 % of the path
is mis-assigned.  The namespace for each G4 name is found
by scanning the release's geometry XMLs rather than guessed.

WHAT IS NOT COVERED
-------------------
r < TrackBeamR1 = 3.15 cm.  The beam pipe lives in
Geometry/CMSCommonData/data/beampipe.xml, which SURVIVES the tracker cull in
runToyGeomCheck.py, so those volumes are already present and must not be
duplicated here.  The traversal below 3.15 cm (Vacuum, Beryllium, Air) is
reported and then dropped.

Fixed z structure.  Coaxial cylinders are faithful at eta = 0.30, where
everything crossed is barrel.  At larger |eta| the real material includes
fixed-z disks a cylinder cannot represent, so the construction is per-eta.

THE FIXED POINT.  Building the toy from the reference's traversal and then
re-running the reference through the toy does not give the same trajectory: the
material moved slightly, so the energy loss and hence the bend differ, so it
crosses slightly different radii.  The radii are taken from the real traversal
and the planes are then derived from the helix ANALYTICALLY (as
gen_toy_config.py does), and the result is VERIFIED rather than assumed -- see
`--verify`, which compares the toy's own dedxdbg dump against the real one.

usage:
    python gen_toy_realmat.py dbg_real.log                       # table only
    python gen_toy_realmat.py dbg_real.log --write               # + xml + planes
    python gen_toy_realmat.py dbg_real.log --verify dbg_toy.log  # closure of the loop
"""
import argparse
import glob
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")

# tracker envelope: Geometry/CMSCommonData/data/cms.xml, TrackBeamR1 = 3.15*cm.
# RMAX/DZ match gen_toy_config.py so the two toys share an envelope.
RMIN, RMAX, DZ = 3.15, 123.0, 194.0

DBG = re.compile(r"\[dedxdbg\] iter=(\S+) len=(\S+) dEdx=(\S+) mat=(\S+) "
                 r"rho=(\S+) r=(\S+) accum=(\S+)")


def parse(path):
    """[(len, dedx, mat, rho, r_end, accum)] in cm / MeV-cm units, in order."""
    out = []
    for line in open(path):
        m = DBG.search(line)
        if m:
            out.append((float(m.group(2)), float(m.group(3)), m.group(4),
                        float(m.group(5)), float(m.group(6)), float(m.group(7))))
    if not out:
        raise SystemExit(f"no [dedxdbg] lines in {path} -- was CVH_DEDX_DEBUG=1 set?")
    return out


def volumes(steps):
    """Coalesce consecutive same-material steps into (mat, rho, r_in, r_out, path)."""
    vols, rprev = [], 0.0
    for length, _dedx, mat, rho, r, _accum in steps:
        if vols and vols[-1]["mat"] == mat:
            vols[-1]["r_out"] = r
            vols[-1]["path"] += length
            vols[-1]["nstep"] += 1
        else:
            vols.append(dict(mat=mat, rho=rho, r_in=rprev, r_out=r,
                             path=length, nstep=1))
        rprev = r
    return vols


def clip(vols, rmin=RMIN, rmax=RMAX):
    """Keep only what the toy is responsible for: rmin <= r <= rmax."""
    kept, dropped = [], []
    for v in vols:
        if v["r_out"] <= rmin or v["r_in"] >= rmax:
            dropped.append(v)
            continue
        w = dict(v)
        if w["r_in"] < rmin:
            # scale the path by the surviving radial fraction; the incidence
            # factor is smooth over a single volume so this is exact to O(dr^2)
            f = (w["r_out"] - rmin) / (w["r_out"] - w["r_in"])
            w["path"] *= f
            w["r_in"] = rmin
            w["clipped"] = True
        if w["r_out"] > rmax:
            f = (rmax - w["r_in"]) / (w["r_out"] - w["r_in"])
            w["path"] *= f
            w["r_out"] = rmax
            w["clipped"] = True
        kept.append(w)
    return kept, dropped


def material_namespaces():
    """G4 material name -> DDD namespace, by scanning the release's XMLs.

    A G4 material carries only the bare name, but `<rMaterial>` needs
    `namespace:name`, and the namespace is the MaterialSection label.  Scanned
    rather than hard-coded so a release change shows up as a missing key here
    instead of a DDD parse error inside cmsRun.
    """
    base = os.environ.get("CMSSW_RELEASE_BASE", "")
    if not base:
        raise SystemExit("CMSSW_RELEASE_BASE unset -- run cmsenv first")
    pat = re.compile(r'<(?:ElementaryMaterial|CompositeMaterial)\s+name="([^"]+)"')
    sec = re.compile(r'<MaterialSection\s+label="([^"]+)"')
    # `materials.xml` first so that a name defined in several files resolves to
    # the common one; trackermaterial.xml next because it is explicitly KEPT
    # through the tracker cull in runToyGeomCheck.py.
    order = ("CMSCommonData/data/materials.xml",
             "TrackerCommonData/data/trackermaterial.xml")
    files = [os.path.join(base, "src/Geometry", o) for o in order]
    files += sorted(glob.glob(os.path.join(base, "src/Geometry/*/data/**/*.xml"),
                              recursive=True))
    ns = {}
    for f in files:
        try:
            s = open(f, errors="ignore").read()
        except OSError:
            continue
        m = sec.search(s)
        if not m:
            continue
        label = m.group(1)[:-4] if m.group(1).endswith(".xml") else m.group(1)
        for name in pat.findall(s):
            ns.setdefault(name, label)
    return ns


def helix_frames(radii, pt, eta, phi0, B, q):
    """Plane origin / normal / u-axis at each radius, from the reference helix.

    Identical construction to gen_toy_config.py: the frame is a property of the
    SURFACE, fixed by the reference, and must be the same object the model side
    is handed.
    """
    R = pt / (0.3 * B) * 100.0
    org, nrm, uu = [], [], []
    for r in radii:
        a = math.asin(min(r / (2 * R), 1.0))
        phi = phi0 - q * a
        z = 2 * R * a * math.sinh(eta)
        org += [r * math.cos(phi), r * math.sin(phi), z]
        nrm += [math.cos(phi), math.sin(phi), 0.0]
        uu += [-math.sin(phi), math.cos(phi), 0.0]
    return org, nrm, uu


def emit_xml(vols, ns, path, split_at=()):
    """One cylinder per volume, with the scored volumes SPLIT at their plane.

    WHY THE SPLIT IS NOT OPTIONAL. ToyStateNtuplizer records the state at the
    post-step point of a step that ended on a GEOMETRY BOUNDARY -- that is what
    makes it exact, with no interpolation and no sagitta error. A plane at the
    mid-plane of a sensor is not a boundary, so Geant4 never terminates a step
    there and the watcher never fires: the first attempt at this geometry ran to
    completion, wrote a well-formed ntuple, and recorded ZERO crossings.

    So each scored volume is emitted as TWO shells meeting at the plane, of
    IDENTICAL material -- pure bookkeeping, no physics, exactly the trick the
    layered toy uses for its own scoring surfaces. It also makes the toy MORE
    like the real geometry, not less: there the propagation terminates on each
    target surface too, so the real sensor is split at the same place.
    """
    edges = []
    for v in vols:
        cut = [r for r in split_at if v["r_in"] + 1e-9 < r < v["r_out"] - 1e-9]
        prev = v["r_in"]
        for r in sorted(cut):
            edges.append(dict(v, r_in=prev, r_out=r))
            prev = r
        edges.append(dict(v, r_in=prev, r_out=v["r_out"]))
    sol, lp, pos, missing = [], [], [], []
    for i, v in enumerate(edges):
        name = f"RM{i:03d}"
        space = ns.get(v["mat"])
        if space is None:
            missing.append(v["mat"])
            space = "materials"
        sol.append(f'  <Tubs name="{name}" rMin="{v["r_in"]:.5f}*cm" '
                   f'rMax="{v["r_out"]:.5f}*cm" dz="{DZ}*cm" '
                   f'startPhi="0*deg" deltaPhi="360*deg"/>')
        lp.append(f'  <LogicalPart name="{name}" category="unspecified">\n'
                  f'   <rSolid name="tracker:{name}"/>\n'
                  f'   <rMaterial name="{space}:{v["mat"]}"/>\n'
                  f'  </LogicalPart>')
        pos.append('  <PosPart copyNumber="1">\n'
                   '   <rParent name="tracker:Tracker"/>\n'
                   f'   <rChild name="tracker:{name}"/>\n  </PosPart>')
    if missing:
        raise SystemExit("no DDD namespace found for: " + ", ".join(sorted(set(missing))))
    hdr = ('<?xml version="1.0"?>\n<DDDefinition xmlns="http://www.cern.ch/cms/DDL" '
           'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
           'xsi:schemaLocation="http://www.cern.ch/cms/DDL '
           '../../../DetectorDescription/Schema/DDLSchema.xsd">\n')
    # NOTE: no ' - - ' inside this comment. XML forbids a double hyphen there and
    # it has cost this project two fatal parse errors.
    body = f"""
 <!-- Toy tracker built from the REAL material sequence along the reference
      trajectory, by Analysis/HitAnalyzer/test/gen_toy_realmat.py. One coaxial
      cylinder per traversed volume, real material by name. Do not hand edit. -->
 <SolidSection label="tracker.xml">
  <Polycone name="Tracker" startPhi="0*deg" deltaPhi="360*deg">
   <ZSection z="-[cms:TrackBeamZ2]" rMin="[cms:TrackBeamR2]" rMax="[cms:TrackCalorR]"/>
   <ZSection z="-[cms:TrackBeamZ1]" rMin="[cms:TrackBeamR1]" rMax="[cms:TrackCalorR]"/>
   <ZSection z="[cms:TrackBeamZ1]"  rMin="[cms:TrackBeamR1]" rMax="[cms:TrackCalorR]"/>
   <ZSection z="[cms:TrackBeamZ2]"  rMin="[cms:TrackBeamR2]" rMax="[cms:TrackCalorR]"/>
  </Polycone>
{chr(10).join(sol)}
 </SolidSection>
 <LogicalPartSection label="tracker.xml">
  <LogicalPart name="Tracker" category="unspecified">
   <rSolid name="tracker:Tracker"/>
   <rMaterial name="materials:Air"/>
  </LogicalPart>
{chr(10).join(lp)}
 </LogicalPartSection>
 <PosPartSection label="tracker.xml">
{chr(10).join(pos)}
 </PosPartSection>
"""
    open(path, "w").write(hdr + body + "\n</DDDefinition>\n")
    return len(edges)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", help="CVH_DEDX_DEBUG=1 output of the REAL-geometry reference")
    ap.add_argument("--pt", type=float, default=3.0)
    ap.add_argument("--eta", type=float, default=0.30)
    ap.add_argument("--phi", type=float, default=0.70)
    ap.add_argument("--bfield", type=float, default=3.8)
    ap.add_argument("--charge", type=float, default=-1.0)
    ap.add_argument("--score", default="Silicon",
                    help="material whose volumes carry a scoring plane")
    ap.add_argument("--plane-at", choices=("mid", "entry"), default="mid",
                    help="mid = the sensor MID-plane, production's own target "
                         "surface (surfacemapD_[detid]), which is NOT a material "
                         "boundary; entry = the inner face, what the layered toy "
                         "and the real clean-prop study both use. The choice is "
                         "not cosmetic: a target on a boundary makes the "
                         "propagation end in a degenerate sliver that Geant4 "
                         "attributes to the FAR side, which is the whole reason "
                         "dEdxlast needed a step-length floor.")
    ap.add_argument("--final-plane", action="store_true",
                    help="append a scoring plane at the OUTER EDGE of the last "
                         "kept volume. The traversal has one fewer Silicon "
                         "volume than the real geometry has targets, because the "
                         "propagation terminates ON each target surface and so "
                         "never enters the last sensor -- which leaves the real "
                         "geometry's outermost LEG with no counterpart here. That "
                         "leg is the most anomalous one in the per-leg "
                         "decomposition, so it needs a "
                         "control. The edge is already a geometry boundary, so no "
                         "split is required and the watcher fires there.")
    ap.add_argument("--write", action="store_true", help="write the xml and the plane file")
    # DDD takes the namespace from the FILE NAME, and the toy has to define
    # tracker:Tracker (cmsTracker.xml positions it into cms:CMSE), so the
    # basename must be tracker.xml. Hence its own directory rather than a
    # different name -- and the layered toy at data/tracker.xml is untouched.
    ap.add_argument("--xml", default=os.path.join(DATA, "realmat", "tracker.xml"))
    ap.add_argument("--planes", default="",
                    help="default: toyPlanes_<geometry dir>_pt3.py next to "
                         "this script, which is the name the runners derive "
                         "from the geometry path.")
    ap.add_argument("--verify", default="",
                    help="dedxdbg log of the TOY reference; compare the two traversals")
    args = ap.parse_args()

    steps = parse(args.log)
    vols = volumes(steps)
    kept, dropped = clip(vols)

    print(f"{len(steps)} steps -> {len(vols)} volumes; "
          f"{len(kept)} inside [{RMIN}, {RMAX}] cm, {len(dropped)} outside\n")
    print("dropped (already in the surviving standard geometry, e.g. the beam pipe):")
    for v in dropped:
        print(f"    {v['mat']:<24} r {v['r_in']:8.4f} - {v['r_out']:8.4f} cm  "
              f"rho {v['rho']:9.5f}")

    print(f"\n{'#':>3} {'material':<24} {'r_in':>9} {'r_out':>9} {'dr':>8} "
          f"{'path':>8} {'rho':>9} {'rho*d':>9} {'steps':>6}")
    xtot = 0.0
    for i, v in enumerate(kept):
        dr = v["r_out"] - v["r_in"]
        x = v["rho"] * v["path"]
        xtot += x
        print(f"{i:3d} {v['mat']:<24} {v['r_in']:9.4f} {v['r_out']:9.4f} "
              f"{dr:8.4f} {v['path']:8.4f} {v['rho']:9.5f} {x:9.5f} {v['nstep']:6d}")
    # PARSE CHECK. path/dr must be the incidence factor 1/(sin(theta) cos(alpha))
    # with sin(alpha) = r/2R -- a smooth function of r and nothing else. If the
    # coalescing were mis-pairing r_in with the wrong step this would scatter.
    R = args.pt / (0.3 * args.bfield) * 100.0
    sth = 1.0 / math.cosh(args.eta)
    worst = 0.0
    for v in kept:
        dr = v["r_out"] - v["r_in"]
        if dr <= 0:
            continue
        rm = 0.5 * (v["r_in"] + v["r_out"])
        pred = 1.0 / (sth * math.cos(math.asin(min(rm / (2 * R), 1.0))))
        worst = max(worst, abs(v["path"] / dr / pred - 1.0))
    print(f"\n    parse check, path/dr against 1/(sin(theta) cos(alpha)): "
          f"worst deviation {worst:.2e}")
    print(f"    total mass thickness along the path: {xtot:.4f} g/cm2")
    print(f"    (the layered toy is 15 x 1 mm x 9 g/cm3 = 13.50 g/cm2 radially)")

    scored = [v for v in kept if v["mat"] == args.score]
    radii = ([0.5 * (v["r_in"] + v["r_out"]) for v in scored] if args.plane_at == "mid"
             else [v["r_in"] for v in scored])
    if args.final_plane and kept:
        redge = kept[-1]["r_out"]
        if not radii or redge > radii[-1] + 1e-6:
            radii.append(redge)
            print(f"\n    + final plane at the traversal's outer edge, "
                  f"r = {redge:.4f} cm (material before it: {kept[-1]['mat']})")
    where = ("MID-plane (production's own target surface)" if args.plane_at == "mid"
             else "INNER FACE (what the layered toy and the real study use)")
    print(f"\n{len(scored)} scoring planes at the {where} of each "
          f"'{args.score}' volume:")
    print("    " + ", ".join(f"{r:.4f}" for r in radii))
    print("\nwatcher radii line for the runner:")
    print("    radii=cms.vdouble(" + ", ".join(f"{r:.4f}" for r in radii) + "),")

    if args.verify:
        vt, _ = clip(volumes(parse(args.verify)))
        n = min(len(kept), len(vt))
        # TWO differences are EXPECTED and are not failures.
        #
        # (1) The toy's last volume is the one holding the outermost scoring
        #     plane, and the propagation ENDS on that plane. Under --plane-at mid
        #     it therefore traverses HALF of it. Interior sensors are unaffected:
        #     their two halves are consecutive steps of the same material and
        #     coalesce back into one volume.
        # (2) The real traversal continues past the last scored sensor to the
        #     next target; the toy has no target there, so the real geometry's
        #     trailing volumes have no counterpart.
        last_expected = 0.5 if args.plane_at == "mid" else 1.0
        print(f"\n=== VERIFY: toy traversal against the real one")
        print(f"    volumes  real {len(kept)}   toy {len(vt)}   compared {n}")
        print(f"    {'#':>3} {'material':<26} {'x_real':>9} {'x_toy':>9} {'ratio':>8}")
        worst, bad = 0.0, 0
        for i in range(n):
            a, b = kept[i], vt[i]
            xa, xb = a["rho"] * a["path"], b["rho"] * b["path"]
            rat = xb / xa if xa else float("nan")
            islast = (i == n - 1)
            tol = 0.01 if not islast else 0.02
            ref = last_expected if islast else 1.0
            off = a["mat"] != b["mat"] or (xa and abs(rat - ref) > tol)
            if off or islast:
                tag = ("   <-- MISMATCH" if off else
                       f"   <-- expected {ref:.2f}: propagation ENDS on this plane")
                print(f"    {i:3d} {a['mat'] + '/' + b['mat']:<26} {xa:9.5f} "
                      f"{xb:9.5f} {rat:8.4f}{tag}")
            if off:
                bad += 1
            if xa and not islast:
                worst = max(worst, abs(rat - 1.0))
        print(f"    worst |x_toy/x_real - 1| over the {n - 1} fully traversed "
              f"volumes: {worst:.2e}")
        print(f"    mismatches: {bad}")
        if len(kept) > n:
            tail = sum(v["rho"] * v["path"] for v in kept[n:])
            print(f"    real-only tail beyond the last scored sensor: "
                  f"{len(kept) - n} volumes, {tail:.5f} g/cm2 "
                  f"({', '.join(v['mat'] for v in kept[n:])})")

    if args.write:
        if not args.planes:
            gdir = os.path.basename(os.path.dirname(os.path.abspath(args.xml)))
            args.planes = os.path.join(HERE, f"toyPlanes_{gdir}_pt3.py")
        ns = material_namespaces()
        nsh = emit_xml(kept, ns, args.xml, split_at=radii)
        print(f"    {len(kept)} volumes -> {nsh} shells "
              f"({nsh - len(kept)} split at a scoring plane; a plane on a volume "
              f"FACE is already a boundary and needs no split)")
        org, nrm, uu = helix_frames(radii, args.pt, args.eta, args.phi,
                                    args.bfield, args.charge)
        f = lambda v: ", ".join("%.6f" % x for x in v)
        open(args.planes, "w").write(
            f"# generated by gen_toy_realmat.py from {os.path.basename(args.log)}; "
            f"do not hand edit\n"
            f"radii = [{f(radii)}]\norigin = [{f(org)}]\n"
            f"normal = [{f(nrm)}]\nuaxis  = [{f(uu)}]\n")
        print(f"\nwrote {args.xml}\nwrote {args.planes}")
        used = sorted({v["mat"] for v in kept})
        print("materials, resolved to DDD namespaces:")
        for m in used:
            print(f"    {ns[m]}:{m}")


if __name__ == "__main__":
    main()
