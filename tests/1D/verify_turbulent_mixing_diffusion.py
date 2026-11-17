#!/usr/bin/env python3
"""
verify_turbulent_mixing_diffusion.py: fvModel turbulentMixingDiffusion
(subgrid products/air mixing in the species and energy of a fluid phase).

Case tests/1D/turbulentMixingDiffusion (blastEulerFoam, one gas phase, closed
planar tube at rest, 1 atm, 200 cells over 0.5 m): pure N2 ("air", YAir = 1)
for x < 0, CO ("products", YProducts = 0) for x > 0; for x > 0.2 m the
products are half O2 (mixture fraction 1, outside the mixing layer).
D(t) = k t from a table Function1 (k = 1250 m^2/s^2), so the diffusion
length 2 sqrt(k t^2/2) is 50 mm (20 cells) at 1 ms.

Runs (copies of the case under tests/1D/turbulent_mixing_diffusion/):
  A  298 K both sides, energy no
  D  as A with k = 0 (no diffusion; reference for the mask check)
  B  298 K | 900 K, energy yes
  C  298 K | 900 K, energy no

Pass criteria:
  A: Y_CO within 1e-3 of 0.5 erfc(-x/(2 sqrt(k t^2/2))) at 0.5 and 1 ms
     (x < 0.15 m); CO and total mass, and total energy, conserved to 1e-9;
     the pressure uniform within 1 Pa.
  mask: the O2 step at x = 0.2 m (z = 1, outside the layer) changes as in D
     (no diffusion) within 1e-4.
  B/C: total energy conserved to 1e-9; with energy diffusion (unity Lewis
     number) the specific energy stays linear in Y_CO within 8% of the
     energy jump (the rest is pdV work of the mixing flow); without it, it
     departs by more than 30%.
"""
import os
import re
import shutil
import subprocess
import sys
from math import erfc, sqrt

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CASE = os.path.join(HERE, "turbulentMixingDiffusion")
WORK = os.path.join(HERE, "turbulent_mixing_diffusion")
K = 1250.0
N = 200
L = 0.25
DX = 2*L/N
X = -L + (np.arange(N) + 0.5)*DX
TIMES = ("0.00025", "0.0005", "0.00075", "0.001")

RUNS = {
    "A": [],
    "D": ["--k", "0"],
    "B": ["--TRight", "900", "--energy", "yes"],
    "C": ["--TRight", "900", "--energy", "no"],
}


def env(cmd):
    return (
        "source $WM_PROJECT_DIR/etc/bashrc >/dev/null 2>&1 && "
        "source $BLAST_DIR/etc/bashrc >/dev/null 2>&1 && " + cmd
    )


def run(name, opts):
    d = os.path.join(WORK, name)
    if os.path.isdir(d):
        shutil.rmtree(d)
    shutil.copytree(CASE, d, ignore=shutil.ignore_patterns(
        "0", "0.*", "log.*", "polyMesh", "postProcessing"))
    return subprocess.Popen(["bash", "-c", env("./Allrun " + " ".join(opts))], cwd=d,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def field(name, t, f):
    s = open(os.path.join(WORK, name, t, f)).read()
    m = re.search(r"nonuniform\s+List<(\w+)>\s*(\d+)\s*\(", s)
    if not m:
        v = re.search(r"internalField\s+uniform\s+([^;]+);", s).group(1)
        return np.zeros(N) if "(" in v else float(v)*np.ones(N)
    k = int(m.group(2))
    if m.group(1) == "vector":
        return np.array([float(v) for v in
                         re.findall(r"\(([^ ()]+) [^()]*\)", s[m.end():])[:k]])
    return np.array(s[m.end():].split(")")[0].split(), dtype=float)[:k]


def main():
    os.makedirs(WORK, exist_ok=True)
    procs = {n: run(n, o) for n, o in RUNS.items()}
    ok = True
    for n, p in procs.items():
        p.wait()
        log = open(os.path.join(WORK, n, "log.blastEulerFoam")).read()
        if "FOAM FATAL" in log or "\nEnd" not in log:
            print(f"run {n}: solver did not finish")
            ok = False
    if not ok:
        print("FAIL: solver runs")
        return 1

    def check(cond, msg):
        nonlocal ok
        print(("  ok   " if cond else "  FAIL ") + msg)
        ok = ok and cond

    # A: erfc solution, conservation, uniform pressure
    inner = X < 0.15
    for t in ("0.0005", "0.001"):
        S = K*float(t)**2/2
        ex = np.array([0.5*erfc(-x/(2*sqrt(S))) for x in X])
        err = np.abs(field("A", t, "CO.gas") - ex)[inner].max()
        check(err < 1e-3, f"A t = {t}: max |Y_CO - erfc| = {err:.2e} (< 1e-3)")
    for n in ("A", "B", "C"):
        ar = [field(n, t, "alphaRho.gas") for t in TIMES]
        co = [np.sum(a*field(n, t, "CO.gas")) for a, t in zip(ar, TIMES)]
        E = [np.sum(field(n, t, "alphaRhoE.gas")) for t in TIMES]
        dm = max(abs(np.sum(a)/np.sum(ar[0]) - 1) for a in ar)
        dco = max(abs(c/co[0] - 1) for c in co)
        dE = max(abs(e/E[0] - 1) for e in E)
        check(max(dm, dco, dE) < 1e-9,
              f"{n}: mass, CO mass, energy conserved to {dm:.1e}, {dco:.1e}, {dE:.1e} (< 1e-9)")
    p = field("A", "0.001", "p")
    check(p.max() - p.min() < 1.0, f"A: pressure range {p.max() - p.min():.2f} Pa (< 1)")

    # Mask: O2 step outside the mixing layer, against no diffusion
    dO2 = np.abs(field("A", "0.001", "O2.gas") - field("D", "0.001", "O2.gas")).max()
    check(dO2 < 1e-4, f"mask: O2 at the z = 1 step differs from D = 0 by {dO2:.1e} (< 1e-4)")

    # B/C: energy follows the species with unity Lewis number
    dev = {}
    for n in ("B", "C"):
        U = field(n, "0.001", "U.gas")
        e = field(n, "0.001", "alphaRhoE.gas")/field(n, "0.001", "alphaRho.gas") - 0.5*U**2
        Y = field(n, "0.001", "CO.gas")
        eL, eR = e[0], e[int(0.75*N)]  # pure N2 and pure CO (x = 0.125 m)
        dev[n] = np.abs(e - (eL + (eR - eL)*Y))[inner].max()/(eR - eL)
    check(dev["B"] < 0.08, f"B (energy yes): max |e - linear in Y_CO| = {dev['B']*100:.1f}% of the jump (< 8%)")
    check(dev["C"] > 0.30, f"C (energy no): {dev['C']*100:.1f}% of the jump (> 30%, energy not diffused)")

    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
