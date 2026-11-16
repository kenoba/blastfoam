#!/usr/bin/env python3
"""Afterburning with air by O2 consumption: O2 used since t = 0 and the
fraction of the products fuel burnt, f = O2 used/(O2 used + O2 still needed to
oxidise all C and H to CO2 and H2O). The domain O2 is from alpha.gas*rho.gas
*O2.gas*Vc (Vc written; time 0 from the initial fields).
usage: o2burn.py <case> [<time> ...]   (default: all written times)"""
import os
import re
import sys

import numpy as np

case = sys.argv[1]


def field(t, name, n=None):
    s = open(os.path.join(case, t, name)).read()
    m = re.search(r"internalField\s+uniform\s+([^;]+);", s)
    if m:
        return float(m.group(1)) * np.ones(n)
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    k = int(m.group(1))
    return np.array(s[m.end():].split(")")[0].split(), dtype=float)[:k]


# O2 needed per mol of species (to CO2, H2O) and molar masses [g/mol]
need = {"CO": 0.5, "C(gr)": 1.0, "H2": 0.5, "CH4": 2.0, "C2H6": 3.5, "H": 0.25,
        "OH": -0.25, "CH3": 1.75, "CH2O": 1.0, "HCO": 0.75, "CH3O": 1.25,
        "O": -0.5, "HO2": -0.25, "H2O2": -0.5, "NH3": 0.75}
W = {"CO": 28.010, "C(gr)": 12.011, "H2": 2.016, "CH4": 16.043, "C2H6": 30.069,
     "H": 1.008, "OH": 17.007, "CH3": 15.035, "CH2O": 30.026, "HCO": 29.018,
     "CH3O": 31.034, "O": 15.999, "HO2": 33.007, "H2O2": 34.015, "NH3": 17.031,
     "O2": 31.998}

times = sys.argv[2:] or sorted(
    [d for d in os.listdir(case) if re.fullmatch(r"[0-9.e-]+", d) and d != "0"],
    key=float)

V0 = field("0", "Vc")
n = len(V0)
O20 = float(np.sum(field("0", "alpha.gas", n) * field("0", "rho.gas", n)
                   * field("0", "O2.gas", n) * V0))
print(f"{'t [ms]':>7} {'O2 used':>9} {'O2 needed':>10} {'burnt':>7}   (mol, domain)")
for t in times:
    if not os.path.isfile(os.path.join(case, t, "Vc")):
        continue
    V = field(t, "Vc")
    ar = field(t, "alpha.gas", len(V)) * field(t, "rho.gas", len(V))
    mol = lambda s: float(np.sum(ar * field(t, s + ".gas", len(V)) * V)) / W[s] * 1e3
    used = O20 / W["O2"] * 1e3 - mol("O2")
    req = sum(need[s] * mol(s) for s in need
              if os.path.isfile(os.path.join(case, t, s + ".gas")))
    print(f"{float(t) * 1e3:7.3f} {used:9.3e} {req:10.3e} {used / (used + req) * 100:6.2f}%")
