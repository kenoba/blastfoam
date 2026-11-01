#!/usr/bin/env python3
"""Products/air contact diagnostics of the 2D axisymmetric cases (G3b).

For each time directory (reconstructed, or processor*/<time> when not
reconstructed) with the cell centres C and volumes (V or Vc; write them
with postProcess -func writeCellCentres / writeCellVolumes):

- products fraction of each cell, f = alpha.compBAl + alpha.gas*phi/phi_p,
  where phi = C/(C + N) is the carbon/(carbon + nitrogen) element mass ratio
  of the gas (phi_p = 0.449 in the products, 0 in air; conserved by the
  chemistry); f is normalised by 1 - alpha.Al;
- on rays at angles psi from the symmetry plane (psiMin to psiMax), the
  contact radius R(psi) = r_core + int_{r > r_core} f dr, r_core the outer
  edge of the last cell with f > 0.99, and the mixing width
  W(psi) = int 4 f (1 - f) dr;
- the amplitude of each mode n over the window,
  a_n = 2/(psiMax - psiMin) int (R - <R>) cos(n psi), and the ray-averaged
  R and W;
- the gas CO and H2 mass in the whole domain.

Usage: contact_modes.py <case> [--times t1 t2 ...] [--modes 24 48]
                        [--psimin 22.5] [--psimax 67.5] [--bin 0.5]
                        [--out file.csv] [--profile file.txt]
The window (default 22.5-67.5 deg) should span a multiple of the
half-period 180/n deg of the modes, and stays clear of the axis (90 deg)
and the symmetry plane (0 deg), where the fireball is mesh-dependent
(axis study, G2b).
"""
import argparse
import glob
import os
import re

import numpy as np

WEDGE = np.radians(5.0)
PHI_P = 0.449
# Atoms of C and N and molar mass [kg/kmol] of the gas species that hold them
CN = {"CO": (1, 0, 28.010), "CO2": (1, 0, 44.009), "CH4": (1, 0, 16.043),
      "C2H6": (2, 0, 30.069), "C(gr)": (1, 0, 12.011), "CH2O": (1, 0, 30.026),
      "CH3": (1, 0, 15.035), "CH3O": (1, 0, 31.034), "HCO": (1, 0, 29.018),
      "N2": (0, 2, 28.013), "NO": (0, 1, 30.006), "N": (0, 1, 14.007),
      "NH3": (0, 1, 17.031)}
HDR = re.compile(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*(\d+)\s*\(")
UNI = re.compile(r"internalField\s+uniform\s+([^;]+);")


def read(path, n=None):
    text = open(path).read()
    m = HDR.search(text)
    if m:
        k = int(m.group(2))
        body = text[m.end():text.index("\n)", m.end())]
        if m.group(1) == "scalar":
            return np.array(body.split(), float)
        return np.array(re.sub(r"[()]", " ", body).split(), float).reshape(k, 3)
    u = UNI.search(text).group(1).strip()
    if u.startswith("("):
        return np.tile(np.array(u.strip("()").split(), float), (n, 1))
    return np.full(n, float(u))


def load(dirs):
    """Concatenate the fields of the given time directories."""
    out = {}
    for d in dirs:
        C = read(os.path.join(d, "C"))
        n = len(C)
        Vf = os.path.join(d, "V") if os.path.exists(os.path.join(d, "V")) \
            else os.path.join(d, "Vc")
        f = {"C": C, "V": read(Vf, n)}
        for name in ("alpha.compBAl", "alpha.gas", "alpha.Al", "alphaRho.gas"):
            f[name] = read(os.path.join(d, name), n)
        for s in CN:
            p = os.path.join(d, s + ".gas")
            f[s] = read(p, n) if os.path.exists(p) else np.zeros(n)
        f["H2"] = read(os.path.join(d, "H2.gas"), n)
        for k, v in f.items():
            out.setdefault(k, []).append(v)
    return {k: np.concatenate(v) for k, v in out.items()}


def analyse(f, modes, psimin, psimax, dbin):
    C, V = f["C"], f["V"]
    r = np.linalg.norm(C, axis=1)
    psi = np.degrees(np.arctan2(C[:, 1], np.hypot(C[:, 0], C[:, 2])))
    wc = sum(f[s]*c*12.011/M for s, (c, nn, M) in CN.items())
    wn = sum(f[s]*nn*14.007/M for s, (c, nn, M) in CN.items())
    phi = np.where(wc + wn > 1e-12, wc/np.maximum(wc + wn, 1e-30), 0.0)
    fp = (f["alpha.compBAl"] + f["alpha.gas"]*np.clip(phi/PHI_P, 0, 1)) \
        / np.maximum(1 - f["alpha.Al"], 1e-6)
    fp = np.clip(fp, 0, 1)
    # Rays every dbin degrees: the cells a ray crosses are those whose
    # centre lies within half a cell size of it (quad cells of side
    # sqrt(area), area = V/(distance from the axis*wedge angle)); along the
    # ray each contributes its side as path length.
    rho_ax = np.hypot(C[:, 0], C[:, 2])
    side = np.sqrt(V/(np.maximum(rho_ax, 1e-9)*WEDGE))
    rmax = 1.1*r[fp > 0.01].max()
    mid = np.arange(psimin + 0.5*dbin, psimax, dbin)
    R, W = np.zeros(len(mid)), np.zeros(len(mid))
    for i, p in enumerate(mid):
        dist = r*np.abs(np.sin(np.radians(psi - p)))
        s = (dist < 0.5*side) & (np.abs(psi - p) < 45) & (r < rmax)
        o = np.argsort(r[s])
        rr, ff, dd = r[s][o], fp[s][o], side[s][o]
        # Contact: the last cell with f > 0.99 from the inside plus the
        # products path length beyond it
        core = np.where(ff > 0.99)[0]
        k = core[-1] if len(core) else 0
        R[i] = rr[k] + 0.5*dd[k] + np.sum(ff[k + 1:]*dd[k + 1:])
        W[i] = np.sum(4*ff*(1 - ff)*dd)
    th = np.radians(mid)
    span = np.radians(psimax - psimin)
    amps = {n: 2/span*np.sum((R - R.mean())*np.cos(n*th))*np.radians(dbin)
            for n in modes}
    full = 2*np.pi/WEDGE
    mCO = full*np.sum(f["alphaRho.gas"]*f["CO"]*V)
    mH2 = full*np.sum(f["alphaRho.gas"]*f["H2"]*V)
    return mid, R, W, amps, mCO, mH2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case")
    ap.add_argument("--times", nargs="*")
    ap.add_argument("--modes", nargs="*", type=int, default=[24, 48])
    ap.add_argument("--psimin", type=float, default=22.5)
    ap.add_argument("--psimax", type=float, default=67.5)
    ap.add_argument("--bin", type=float, default=0.5)
    ap.add_argument("--out")
    ap.add_argument("--profile", help="write R(psi), W(psi) of each time here")
    a = ap.parse_args()

    def is_time(name):
        try:
            float(name)
            return True
        except ValueError:
            return False

    times = a.times
    if not times:
        names = os.listdir(a.case)
        if not any(is_time(n) and n != "0" for n in names):
            names = os.listdir(os.path.join(a.case, "processor0"))
        times = sorted({n for n in names if is_time(n) and n != "0"}, key=float)
    rows = []
    for t in times:
        d = os.path.join(a.case, t)
        dirs = [d] if os.path.exists(os.path.join(d, "C")) else \
            sorted(glob.glob(os.path.join(a.case, "processor*", t)))
        if not dirs or not all(os.path.exists(os.path.join(x, "C")) for x in dirs):
            print(f"{t}: no cell centres, skipped")
            continue
        mid, R, W, amps, mCO, mH2 = analyse(load(dirs), a.modes, a.psimin, a.psimax, a.bin)
        rows.append((float(t), R.mean(), R.std(), W.mean(), mCO, mH2,
                     *[amps[n] for n in a.modes]))
        print(f"t={float(t)*1e6:7.1f} us  <R>={R.mean()*1e3:7.2f} mm  "
              f"sd(R)={R.std()*1e3:6.3f} mm  <W>={W.mean()*1e3:6.2f} mm  "
              f"CO={mCO*1e3:7.2f} g  H2={mH2*1e3:6.3f} g  " +
              "  ".join(f"a{n}={amps[n]*1e3:+.3f} mm" for n in a.modes))
        if a.profile:
            with open(a.profile, "a") as fo:
                for p, rr, ww in zip(mid, R, W):
                    fo.write(f"{t} {p} {rr} {ww}\n")
    if a.out:
        with open(a.out, "w") as fo:
            fo.write("t,R_mean,R_sd,W_mean,m_CO,m_H2," +
                     ",".join(f"a{n}" for n in a.modes) + "\n")
            for row in rows:
                fo.write(",".join(f"{v:.6g}" for v in row) + "\n")


if __name__ == "__main__":
    main()
