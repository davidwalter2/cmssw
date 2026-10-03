#!/usr/bin/env python3
"""Write the nuclear-elastic (hadElastic) kernel tables that `cvhcf` reads:
`cvhcf_nucel_v2.bin`, one record per (species, ELEMENT, momentum node).

WHAT THE TABLES ARE.  For each of pi+, pi-, K+, K-, p, pbar and each element of
the job's Geant4 material table, Geant4's own elastic model and cross-section
dataset for that species (G4ElasticHadrNucleusHE/G4BGGPionElasticXS for pions,
G4HadronElastic/Glauber-Gribov for kaons, G4ChipsElasticModel/
G4BGGNucleonElasticXS for protons, G4AntiNuclElastic (G4HadronElastic below
100 MeV)/AntiAGlauber for antiprotons) were sampled with 5e5 collisions per
momentum node by `nucel_g4driver` at the element's own Z and atomic mass, and
reduced to
  * the rate mu = Sigma/rho [cm^2/g] on a fine momentum grid,
  * the projected single-collision CF of the deflection, g(u) = E[J0(u theta)],
  * the distribution of the kinetic energy the projectile gives to the
    recoiling nucleus, dE [MeV],
  * the moments <theta>, <theta^2>, <dE>, <dE^2>,
  * the joint (theta, dE) law of the same collisions on 200 log-theta bins.
Built by `calibration_studies/resolution/ksclosure/nucel/nucel_tables.py`
(build + assemble); that module's `Table` class is the offline reference
implementation of every evaluation rule below and reads THIS file, so the
references -- the fit-level `ks_nucel_cf.step_family` and the
clean-propagation `cf_nucel_exact.rows_exponent` -- and the in-maker port
evaluate identical numbers.

WHY THIS FORM.  The in-maker family needs, per MS step s of block b, the
angular kernel at the argument w_b * tau * (momentum rescale) -- a continuous,
block-dependent multiple of the 64-point tau grid (CvhCfExponents.h) -- and
the recoil kernel at wq_s * tau with a per-step ionisation weight.  A table on
the tau grid itself cannot serve a continuous weight.  Of the two remaining
options:
  * angular: the CF g(u) on a log-u grid.  It is smooth in ln u; a
    Catmull-Rom interpolation of the 30-per-decade nodes reproduces the
    sampled CF on a 2.3x denser grid to 4.9e-4 (max over all 18426 records,
    u <theta> <= 20; median 6e-5), below the ~1e-3 Monte-Carlo noise of the
    sampled CF, and evaluating it is a table lookup.  Tabulating the theta DENSITY instead
    would put a quadrature of J0 (no closed form per bin) over ~1e4 bins x 64
    tau x every step into the maker, and make the answer depend on the
    toolchain's J0 (the reason the Moliere shape tables are shipped as data).
  * recoil: the dE DISTRIBUTION on log bins (partial moments), with the CF
    evaluated per bin in closed form (sin/cos only).  The recoil CF oscillates
    in v wherever the spectrum has structure on the scale 1/v (hydrogen
    targets: a kinematic end point), so a CF tabulated on a v grid is wrong
    between its nodes at the per-cent level; a bin-uniform density is wrong
    by up to 4e-2 at 16 bins per decade; matching each bin's mean and
    variance at 32 bins per decade reproduces the exact sample CF to the
    sample's own noise (max 1.9e-3 over v <dE> <= 20, 48 random buckets), and
    needs no v range.
The compound mixture, the momentum interpolation and both evaluations are
linear in the stored values, so they commute with each other.

FLOAT32.  The per-node records are float32 (the grids, edges, node momenta and
ln mu are float64).  g - 1 is stored rather than g so that the small-u region
keeps its relative precision, and the recoil bins carry moments about the bin
centre so that the in-bin variance is not a float32 difference of two large
numbers.  Measured (`nucel_table_checks.py precision`) on the four gate
compounds, every species and every base node, against the float64 assembly:
max |g32 - g64| = 6.3e-8 over the whole u range; max |h32 - h64| = 7.8e-8 for
v <dE> <= 20 and 3.3e-5 out to v = 1e3/MeV (the phase v * mu_b of recoils of
1e4 MeV) -- against the 7.6e-6 float32 floor of the exported exponents and the
~1e-3 Monte-Carlo noise of the sampled kernels.

BINARY LAYOUT (little-endian, no padding; i4 = int32, f4 = float32, f8 = float64)

  offset 0     char[8]  magic "CVHNUCEL"
               i4       version = 2
               i4       nS     number of species (6)
               i4       nE     number of elements (37)
               i4       nU     angular grid points (241)
               i4       nD     recoil density bins (256)
               i4       nR     rate grid points (561)
               f8[nU]   uGrid     1/rad, uniform in ln u: u_k = u_0 (u_{nU-1}/u_0)^(k/(nU-1))
                                  (0.1 .. 1e7, 30 per decade)
               f8[nD+1] deEdges   MeV, bin edges of the recoil density, uniform in
                                  ln dE (1e-3 .. 1e5, 32 per decade)
               f8[nR]   pRate     GeV/c, rate nodes (0.18 .. 48, uniform in ln p)
               i4[nS]   speciesPdg   211, -211, 321, -321, 2212, -2212
               i4[nE]   elemZ        ascending
               f8[nE]   elemA        g/mole, the atomic mass the element was sampled
                                     with (Geant4's natural-isotope value)
  then, for each species s in order:
               i4       nEdge_s
               i4       nNode_s
               f8[nEdge_s]  edge_s     GeV/c, momenta where the species' model
                                       changes its sampling table (kernels are
                                       discontinuous there; never interpolate across)
               f8[nNode_s]  nodeP_s    GeV/c, ascending kernel nodes
               i4[nNode_s]  nodeSeg_s  segment of each node = number of edges <= nodeP
  then         f8[nS][nE][nR]   lnMu   ln(mu / (cm^2/g)) at pRate
  The records.  RECORD ORDER everywhere below: species s (file order), then
  element e, then node n of species s; nRec = nE * sum_s nNode_s.
  then         f4[nRec][7 + nU]   the fixed part of every record:
       [0] <theta>      rad      } moments of the sampled collisions at nodeP
       [1] <theta^2>    rad^2    }
       [2] <dE>         MeV      }
       [3] <dE^2>       MeV^2    }
       [4] P_low        P(dE < deEdges[0])             } the recoil below the first
       [5] M1_low       E[dE 1{dE < deEdges[0]}] MeV   } bin, carried by its
       [6] M2_low       E[dE^2 1{..}] MeV^2            } partial moments
       [7 .. 7+nU-1]    g(u_k) - 1
  then         i4[nRec][2]   (b0, nb): the recoil bins b0 .. b0+nb-1 are the
                             record's nonzero range (nb = 0: no bin)
  then         f4[sum nb * 3]  per record in order, three blocks of nb values:
                             P_b, D1_b, D2_b = E[(1, dE - c_b, (dE - c_b)^2) 1{deEdges[b] <= dE < deEdges[b+1]}]
                             (MeV powers), c_b = (deEdges[b] + deEdges[b+1]) / 2: partial
                             moments about the bin centre (linear in the distribution, and
                             free of the float32 cancellation of a variance formed from raw
                             moments); bins outside [b0, b0+nb) are zero.
     (P_low + sum_b P_b = 1; no recoil reaches deEdges[nD].)
  then         i4[nRec]      nJ: the record's number of joint bins
  then         f4[sum nJ * 3]  per record in order, three blocks of nJ values:
                             theta_k [rad], dE_k [MeV], w_k -- the k-th nonempty
                             log-theta bin's mean deflection, mean recoil and
                             share of ALL the record's collisions (zero
                             deflections are in no bin); sum_k w_k <= 1.
  then         char[8]  trailer "CVHNUEND"; the file ends there.

EVALUATION (what `nucel_tables.Table` does, and what the port must do)

  Species.  From the maker's configured particle and the track charge, never
  from the mass (p and pbar share a mass, not a model).  Muons: no table.

  Material mixture, per (species, material) and node n.  For the material's
  elements i (mass fraction w_i, atomic mass A_i; elements with w_i = 0
  skipped), r_i(n) = w_i (elemA_i / A_i) exp(lnMu_i interpolated log-log at
  nodeP_n); the record of the material at node n is sum_i r_i K_i(n) /
  sum_i r_i for every record entry K, the recoil bins included (they are
  partial moments, linear in the distribution).  The rate on the rate grid is
  mu_mat(p_r) = sum_i w_i (elemA_i / A_i) exp(lnMu_i(p_r)), stored as its log.
  (The factor elemA_i/A_i: at a fixed per-atom cross section the mass
  attenuation scales as 1/A.)  A material with an element not in the table
  has no tabulated kernel: fail loudly or skip the family for it, and report.

  Rate at momentum p: linear interpolation of ln mu_mat in ln p between the
  two bracketing pRate nodes; constant beyond the ends.  N_s = mu_mat(p_s) *
  xg_s [g/cm^2].

  Node weights at momentum p.  seg = number of edge_s <= p.  Among the nodes
  with nodeSeg == seg (a contiguous range [a, b]): if p <= nodeP_a use node a
  alone; if p >= nodeP_b use node b alone; else the bracketing pair j, j+1
  with f = ln(p/nodeP_j) / ln(nodeP_{j+1}/nodeP_j), weights (1-f, f).

  Angular CF of one node record at argument u >= 0:
     u < u_0:          g = 1 - u^2 <theta^2> / 4
     u > u_{nU-1}:     g = 1 + y[nU-1]
     otherwise:        t = (ln u - ln u_0) / D,  D = (ln u_{nU-1} - ln u_0) / (nU-1)
                       i = clamp(floor(t), 0, nU-2),  s = t - i
                       p0 = y[i-1] (i = 0: 2 y[0] - y[1]),  p1 = y[i],  p2 = y[i+1],
                       p3 = y[i+2] (i = nU-2: 2 y[nU-1] - y[nU-2])
                       g = 1 + 0.5 (2 p1 + (p2 - p0) s + (2 p0 - 5 p1 + 4 p2 - p3) s^2
                                    + (3 p1 - p0 - 3 p2 + p3) s^3)
     with y = the record's g - 1 values.  At momentum p with node weights
     (j, 1-f), (j+1, f):  g(p; x) = (1-f) g_j(x nodeP_j / p) + f g_{j+1}(x nodeP_{j+1} / p)
     (each node rescaled to fixed momentum transfer).

  Recoil CF of one node record at argument v (any sign): each bin is taken
  uniform over mu_b -+ a_b, mu_b = c_b + D1_b/P_b,
  a_b = sqrt(3 max(D2_b/P_b - (D1_b/P_b)^2, 0)) (the bin's own mean and
  variance), so
     h = P_low + i v M1_low - v^2 M2_low / 2
         + sum_{b: P_b > 0} P_b exp(i v mu_b) sinc(v a_b),   sinc(x) = sin(x)/x (1 at 0).
  The family needs h - 1, which is formed without cancellation and
  normalised exactly (P_low + sum_b P_b = 1 drops the constant; the float32
  P's sum to 1 only to their precision):
     h - 1 = i v M1_low - v^2 M2_low / 2 + sum_b P_b [e^{i v mu_b} sinc(v a_b) - 1],
     e^{ix} s - 1 = -2 sin^2(x/2) s + (s - 1) + i sin(x) s,
  with sinc - 1 by its Taylor series below |v a_b| = 1.
  At momentum p: h(p; v) = (1-f) h_j(v) + f h_{j+1}(v), each node's (material)
  record evaluated on its own; no argument rescale: the recoil t/2M is
  momentum independent at fixed momentum transfer.

  Joint bins of a (material) node: the concatenation of its elements' bins,
  element i's weights times its collision share r_i(n) / sum_j r_j(n);
  elements with a share below 1e-5 are left out (trace elements).  At
  momentum p the two nodes' bins with weights (1-f, f), each deflection
  rescaled to fixed momentum transfer (theta p_j / p), the recoil not.  The
  joint term of one collision at angular argument b and recoil argument a is
     J = sum_k W_k [(J0(b theta_k) - 1)(e^{i a X_k} - 1) + (e^{i a X_k} - e^{i a dE_k})],
  X_k = T_eff(dE_k) (the exact 1/p map at the step's E, p) under the exact
  map, dE_k otherwise -- and for a recoil beyond the knock-on channel's p'
  floor (cf_knockon.PMIN_FRAC: the upper node's recoils near the kinematic
  end point can exceed the step's kinetic energy, the recoil being carried at
  fixed momentum transfer).

  The family then is S_ang(tau) = sum_s N_s (g(p_s; w_b tau) - 1),
  S_rec(tau) = sum_s N_s (h(p_s; wq_s tau) - 1) and
  S_jnt(tau) = sum_s N_s J(p_s; w_b tau, wq_s tau) (ksclosure/nucel/ks_nucel_cf.py,
  cf_nucel_exact.rows_exponent).

MODEL TAG.  The format and content id is the file name stem,
"cvhcf_nucel_v2"; `modelTag()` carries "nuc=cvhcf_nucel_v2" when the family
is exported.  Bump the version (and the file name) for any change of
layout, grids, species, elements or sampling.

usage:
    source /work/submit/david_w/ZMass/calibration_studies/setup_env.sh
    python3 make_cvhcf_nucel_tables.py [--npz <assembled table.npz>] [-o cvhcf_nucel_v2.bin]
"""
import argparse
import hashlib
import os
import struct
import sys

import numpy as np

NUCEL = "/work/submit/david_w/ZMass/calibration_studies/resolution/ksclosure/nucel"
NPZ = "/ceph/submit/data/user/d/david_w/ZMass/cvh/nucel_trackfit_260927_v2/nucel_table_v2.npz"
MAGIC = b"CVHNUCEL"
TAIL = b"CVHNUEND"
VERSION = 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", default=NPZ, help="nucel_tables.py assemble output")
    ap.add_argument("-o", "--out", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "cvhcf_nucel_v2.bin"))
    a = ap.parse_args()

    sys.path.insert(0, NUCEL)
    import nucel_tables as nt                                # noqa: E402

    z = np.load(a.npz)
    species = [int(x) for x in z["species"]]
    ug, dee, prate = z["ug"], z["dee"], z["prate"]
    elemZ, elemA = z["elemZ"].astype(np.int32), z["elemA"].astype(np.float64)
    nS, nE, nU, nD, nR = len(species), len(elemZ), len(ug), len(dee) - 1, len(prate)
    assert np.allclose(np.diff(np.log(ug)), np.log(ug[-1] / ug[0]) / (nU - 1), rtol=1e-9, atol=0)
    assert np.all(np.diff(elemZ) > 0)
    recs = []
    with open(a.out, "wb") as f:
        f.write(MAGIC)
        f.write(struct.pack("<6i", VERSION, nS, nE, nU, nD, nR))
        f.write(np.ascontiguousarray(ug, "<f8").tobytes())
        f.write(np.ascontiguousarray(dee, "<f8").tobytes())
        f.write(np.ascontiguousarray(prate, "<f8").tobytes())
        f.write(np.asarray(species, "<i4").tobytes())
        f.write(np.ascontiguousarray(elemZ, "<i4").tobytes())
        f.write(np.ascontiguousarray(elemA, "<f8").tobytes())
        for pdg in species:
            e, p, s = z[f"{pdg}_edges"], z[f"{pdg}_nodes"], z[f"{pdg}_seg"]
            assert np.array_equal(s, np.searchsorted(e, p, side="right"))
            f.write(struct.pack("<2i", len(e), len(p)))
            f.write(np.ascontiguousarray(e, "<f8").tobytes())
            f.write(np.ascontiguousarray(p, "<f8").tobytes())
            f.write(np.ascontiguousarray(s, "<i4").tobytes())
        f.write(np.ascontiguousarray(z["logmu"], "<f8").tobytes())
        # fixed part of every record
        for pdg in species:
            r = np.concatenate([z[f"{pdg}_mom"], z[f"{pdg}_low"], z[f"{pdg}_gm1"]], axis=-1)
            assert r.shape[-1] == 7 + nU
            f.write(np.ascontiguousarray(r, "<f4").tobytes())
        # the recoil bins: nonzero range of every record, then its data
        idx, data = [], []
        for pdg in species:
            pde = nt.centre_moments(z[f"{pdg}_pde"], dee)  # (element, node, 3, nD)
            assert pde.shape[-2:] == (3, nD)
            for e in range(nE):
                for n in range(pde.shape[1]):
                    nz = np.where(pde[e, n, 0] > 0.0)[0]
                    b0, nb = (int(nz[0]), int(nz[-1] - nz[0] + 1)) if len(nz) else (0, 0)
                    idx.append((b0, nb))
                    data.append(pde[e, n, :, b0:b0 + nb].ravel())
            recs.append(np.concatenate([z[f"{pdg}_mom"], z[f"{pdg}_low"], z[f"{pdg}_gm1"],
                                        pde.reshape(pde.shape[:-2] + (-1,))], axis=-1))
        f.write(np.asarray(idx, "<i4").tobytes())
        f.write(np.ascontiguousarray(np.concatenate(data), "<f4").tobytes())
        # the joint bins: count of every record, then its (theta, dE, w)
        jn, jdata, joints = [], [], []
        for pdg in species:
            cnt = z[f"{pdg}_jcnt"]
            th, de, w = z[f"{pdg}_jth"], z[f"{pdg}_jde"], z[f"{pdg}_jw"]
            off = np.concatenate([[0], np.cumsum(cnt.ravel())])
            assert off[-1] == len(th) == len(de) == len(w)
            k = 0
            for e in range(nE):
                for n in range(cnt.shape[1]):
                    a0, a1 = off[k], off[k + 1]
                    jn.append(a1 - a0)
                    jdata.append(np.concatenate([th[a0:a1], de[a0:a1], w[a0:a1]]))
                    joints.append((th[a0:a1], de[a0:a1], w[a0:a1]))
                    k += 1
        f.write(np.asarray(jn, "<i4").tobytes())
        f.write(np.ascontiguousarray(np.concatenate(jdata), "<f4").tobytes())
        f.write(TAIL)

    # read back with the offline reference's reader
    t = nt.read_bin(a.out)
    k = 0
    for s in range(nS):
        assert np.array_equal(t["rec"][s], recs[s].astype(np.float32).astype(np.float64))
        for e in range(nE):
            for n in range(len(t["nodes"][s])):
                for q in range(3):
                    assert np.array_equal(t["joint"][s][e][n][q],
                                          joints[k][q].astype(np.float32).astype(np.float64))
                k += 1
    assert np.array_equal(t["logmu"], z["logmu"]) and np.array_equal(t["ug"], ug)
    n = os.path.getsize(a.out)
    md5 = hashlib.md5(open(a.out, "rb").read()).hexdigest()
    print(f"wrote {a.out}  {n} bytes  md5 {md5}")
    print(f"  species {species}")
    print(f"  elements {nE}: Z {[int(x) for x in elemZ]}")
    print(f"  u grid {nU} [{ug[0]:g}, {ug[-1]:g}]  dE bins {nD} [{dee[0]:g}, {dee[-1]:g}] MeV  "
          f"rate grid {nR} [{prate[0]:g}, {prate[-1]:g}] GeV")
    for pdg in species:
        ed = [round(float(x), 7) for x in z[f"{pdg}_edges"]]
        print(f"  {pdg:6d}: {len(z[f'{pdg}_nodes'])} nodes, edges {ed}")
    print(f"  sampling: {int(z['nsamp'])} collisions per bucket, seed base {int(z['seed'])}")
    for m in z["models"]:
        print(f"  model {m}")


if __name__ == "__main__":
    main()
