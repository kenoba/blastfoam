#!/usr/bin/env python3
"""Rewrite the D(t) table of the turbulentMixingDiffusion fvModel ('mixing')
in <case>/constant/fvModels up to a new end time, with the D(t) of
setup_mix.py (C w dw/dt from the Balakrishnan et al. 2010 mixing-layer width,
decaying as (tau1/tau)^late after tau1). Used to continue a 1D mixing run on
a longer domain (map1d.py) or in 2D (setup_foi2013_2d.py): those cases run
past the endTime of the original table, beyond which the table is held
constant.

usage: extend_mix_table.py <case> --endTime 0.018 [--C 0.005] [--tau1 0.9] [--late 1]
The charge (m_total, R) is read from <case>/charge.txt; C, tau1 and late
default to the values in the existing fvModels comment. The new table is
checked against the old one on its time range."""
import argparse
import os
import re
import sys

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("case")
ap.add_argument("--endTime", type=float, required=True)
ap.add_argument("--C", type=float)
ap.add_argument("--tau1", type=float)
ap.add_argument("--late", type=float)
a = ap.parse_args()

fv = os.path.join(a.case, "constant", "fvModels")
s = open(fv).read()
m = re.search(r"\nmixing\s*\{", s)
if not m:
    sys.exit(f"no 'mixing' fvModel in {fv}")
note = re.search(r"--C ([0-9.eE+-]+) --tau1 ([0-9.eE+-]+) --late ([0-9.eE+-]+)", s)
C = a.C if a.C is not None else float(note.group(1))
tau1 = a.tau1 if a.tau1 is not None else float(note.group(2))
late = a.late if a.late is not None else float(note.group(3))

info = open(os.path.join(a.case, "charge.txt")).read().split()
W = float(info[info.index("m_total") + 1])
r0 = float(info[info.index("R") + 1])
W3 = W ** (1.0 / 3.0)


def width(tau):
    return r0 * max(7.4 * tau, 19.7 * tau ** 1.56)


def D(t):
    """C w dw/dt [m^2/s] at time t [s], as setup_mix.py"""
    if t <= 0:
        return 0.0
    tau = t * 1e3 / W3
    tc = min(tau, tau1)
    h = 1e-4 * tc
    dwdtau = (width(tc + h) - width(tc - h)) / (2 * h) if tc > h else 7.4 * r0
    d = C * width(tc) * dwdtau * 1e3 / W3
    return d * (tau1 / tau) ** late if tau > tau1 else d


block = s[m.start():]
old = np.array(re.findall(r"\(([0-9.eE+-]+) ([0-9.eE+-]+)\)", block), dtype=float)
# compare away from the kink at tau1, where the old table interpolates linearly
far = np.abs(old[:, 0] * 1e3 / W3 - tau1) > 0.02
dev = max((abs(D(t) - d) / d for t, d in old[far] if t > 0), default=0.0)
if dev > 1e-4:
    sys.exit(f"the D formula does not reproduce the existing table (max rel. deviation {dev:.2e})")

ts = sorted(set([0.0] + [a.endTime * (i / 400.0) ** 2 for i in range(1, 401)]))
table = "\n".join(f"            ({t:.6e} {D(t):.6e})" for t in ts)
new_block = re.sub(r"(values\s*\(\n).*?(\n\s*\);)", lambda x: x.group(1) + table + x.group(2),
                   block, count=1, flags=re.S)
s = s[:m.start()] + new_block
s = re.sub(r"// \((D table to [^;]*?; )?setup_mix.py", f"// (D table to {a.endTime:g} s; setup_mix.py", s, count=1)
open(fv, "w").write(s)
print(f"{fv}: C {C:g}, W {W} kg, r0 {r0} m; {len(ts)} rows to {a.endTime:g} s "
      f"(old table reproduced to {dev:.1e}); D(endTime) = {D(a.endTime):.4g} m2/s")
