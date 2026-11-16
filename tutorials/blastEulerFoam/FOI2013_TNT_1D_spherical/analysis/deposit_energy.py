#!/usr/bin/env python3
"""Add energy to the TNT products in a 1D FOI case at one time: a stand-in
for afterburning energy, for a sensitivity of the far-field blast.

    deposit_energy.py <case> <time> <fraction> <charge kg in domain> <e0 J/kg>

The energy fraction*e0*m_charge goes to the gas in proportion to its products
mass, f = Y_C/Y_C,products (the air has no carbon), at fixed density and
composition: T.gas from e_s(T_new) = e_s(T_old) + de (nasa9 sensible energy,
condensed species e = h), T = T.gas where the gas fills the cell, and p, p.gas
scaled with T_new/T_old (Abel-Nobel at fixed rho and composition)."""
import math
import os
import re
import sys

import numpy as np

case, t, frac, mcharge, e0 = sys.argv[1], sys.argv[2], *map(float, sys.argv[3:6])
R = 8.314462618
d = os.path.join(case, t)


def read(name):
    s = open(os.path.join(d, name)).read()
    m = re.search(r"internalField\s+uniform\s+([^;]+);", s)
    if m:
        return s, None, float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    n = int(m.group(1))
    end = s.index(")", m.end())
    return s, (m.start(), s.index(";", end) + 1), np.array(s[m.end():end].split(), dtype=float)[:n]


def write(name, v):
    s, span, _ = read(name)
    if span is None:
        m = re.search(r"internalField\s+uniform\s+[^;]+;", s)
        span = (m.start(), m.end())
    body = "\n".join(f"{x:.12g}" for x in v)
    open(os.path.join(d, name), "w").write(
        s[:span[0]] + f"internalField   nonuniform List<scalar> \n{len(v)}\n(\n{body}\n)\n;" + s[span[1]:])


def field(name, n=None):
    v = read(name)[2]
    return np.full(n, v) if np.ndim(v) == 0 else v


th = open(os.path.join(case, "constant", "thermo.nasa9Gas")).read()
species = re.findall(r'^"?([A-Za-z][\w()]*)"?\s*\n\{', th, re.M)
coef = {}
for sp in species:
    blk = th[th.index("\n" + ("\"%s\"" % sp if "(" in sp else sp) + "\n"):]
    blk = blk[:blk.index("\n}\n")]
    W = float(re.search(r"molWeight\s+([\d.eE+-]+);", blk).group(1))
    Ts = [float(x) for x in re.search(r"Tswitch\s*\(([^)]*)\)", blk).group(1).split()]
    cs = [[float(x) for x in c.split()] for c in
          re.findall(r"\(([-\d.eE+ ]+)\)", blk[blk.index("coeffs"):])]
    cond = re.search(r"condensed\s+(yes|true|on)", blk) is not None
    coef[sp] = (W, Ts, cs, cond)


def es(sp, T):
    """sensible energy [J/kg] relative to 298.15 K"""
    W, Ts, cs, cond = coef[sp]

    def h(T):
        i = min(max(np.searchsorted(Ts, T) - 1, 0), len(cs) - 1)
        a = cs[i]
        return (-a[0]/T**2 + a[1]*math.log(T)/T + a[2] + a[3]*T/2 + a[4]*T**2/3
                + a[5]*T**3/4 + a[6]*T**4/5 + a[7]/T)*R*T/W*1e3
    e = h(T) - h(298.15)
    if not cond:
        e -= R/W*1e3*(T - 298.15)
    return e


n = len(field("Ccx")) if os.path.exists(os.path.join(d, "Ccx")) else None
Tg = field("T.gas")
n = len(Tg)
Y = {sp: field(sp + ".gas", n) for sp in species if os.path.exists(os.path.join(d, sp + ".gas"))}
carbon = {"CO": 1, "C(gr)": 1, "CO2": 1, "CH4": 1, "CH2O": 1, "HCO": 1, "C2H6": 2, "CH3": 1, "CH3O": 1}
YC = sum(Y[s]*carbon[s]*12.011/coef[s][0] for s in carbon if s in Y)
a, rho, V = field("alpha.gas", n), field("rho.gas", n), field("Vc", n)
m = a*rho*V
YCprod = (YC*m).sum()/mcharge
f = np.clip(YC/YCprod, 0, 1)
mprod = (f*m).sum()
de = frac*e0*mcharge/mprod
print(f"charge {mcharge:.6g} kg, products in gas {mprod:.6g} kg, Y_C products {YCprod:.4f}, "
      f"added {frac*e0*mcharge:.6g} J ({frac:.0%} of e0*m), {de/1e6:.3f} MJ/kg of products")


def emix(i, T):
    return sum(Y[s][i]*es(s, T) for s in Y)


Tn = Tg.copy()
for i in np.nonzero(f*de > 1.0)[0]:
    target = emix(i, Tg[i]) + f[i]*de
    T = Tg[i]
    for _ in range(50):
        cv = (emix(i, T + 0.5) - emix(i, T - 0.5))
        dT = (target - emix(i, T))/cv
        T += dT
        if abs(dT) < 1e-6:
            break
    Tn[i] = T
added = sum(m[i]*(emix(i, Tn[i]) - emix(i, Tg[i])) for i in range(n) if Tn[i] != Tg[i])
print(f"check: added {added:.6g} J; T.gas max {Tg.max():.0f} -> {Tn.max():.0f} K, "
      f"products-mass mean {np.sum(f*m*Tg)/mprod:.0f} -> {np.sum(f*m*Tn)/mprod:.0f} K")
r = Tn/Tg
write("T.gas", Tn)
full = a > 0.999
write("T", np.where(full, Tn, field("T", n)*r))
for p in ("p", "p.gas"):
    write(p, field(p, n)*r)
