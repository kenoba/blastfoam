#!/usr/bin/env python3
"""Volume-weighted mapping of a 1D spherical (wedge) blastEulerFoam solution
onto another 1D spherical mesh, possibly longer (cells beyond the source take
the outermost source cell, i.e. ambient air).

Primitive fields only (the solver encodes the conserved fields at start):
- alpha.<phase>, p, rho, T (mixture), scalars without a phase: volume average
- rho.<phase>: alpha-weighted (conserves alpha*rho)
- U.<phase>, T.<phase>, species and other <phase> scalars: mass-weighted
  with alpha.<phase>*rho.<phase>
- U (mixture): mass-weighted with rho
- pTemporalMax: maximum of the overlapping source cells

Usage: map1d.py <source case> <source time> <target case> <target time>
Cell boundaries are taken from the cell centres written by
postProcess -funcs '(writeCellCentres)' (Ccx) in both cases' time dirs."""

import os
import re
import sys

import numpy as np

src, st, tgt, tt = sys.argv[1:5]

HDR = re.compile(r"internalField\s+(nonuniform\s+List<(\w+)>\s*(\d+)\s*\(|uniform\s+)")


def read(path):
    s = open(path).read()
    m = HDR.search(s)
    if m.group(2):                       # nonuniform
        n = int(m.group(3))
        start = m.end()
        if m.group(2) == "scalar":
            # the list ends at the first ")" followed by ";"
            end = re.compile(r"\)\s*;").search(s, start).start()
            vals = np.array(s[start:end].split(), dtype=float)
        else:
            end = re.compile(r"\)\s*\)\s*;").search(s, start).start() + 1
            vals = np.array(re.findall(r"\(([^()]*)\)", s[start:end]), dtype=str)
            vals = np.array([v.split() for v in vals], dtype=float)
        assert len(vals) == n, (path, len(vals), n)
        kind = m.group(2)
        semi = s.index(";", end)
        return s, kind, vals, m.start(), semi + 1
    # uniform
    semi = s.index(";", m.end())
    txt = s[m.end():semi].strip()
    if txt.startswith("("):
        return s, "vector", np.array(txt.strip("()").split(), dtype=float), m.start(), semi + 1
    return s, "scalar", float(txt), m.start(), semi + 1


def values(path, n):
    s, kind, v, a, b = read(path)
    if np.ndim(v) == 0 or (kind == "vector" and np.ndim(v) == 1):
        v = np.tile(v, (n, 1)) if kind == "vector" else np.full(n, v)
    return kind, v


def write(path, template, kind, v):
    s, _, _, a, b = read(template)
    if kind == "scalar":
        body = "\n".join(f"{x:.12g}" for x in v)
    else:
        body = "\n".join("(" + " ".join(f"{c:.12g}" for c in x) + ")" for x in v)
    new = (s[:a] + f"internalField   nonuniform List<{kind}> \n{len(v)}\n(\n"
           + body + "\n)\n;" + s[b:])
    open(path, "w").write(new)


def edges(case, time):
    """Radial cell edges: the distinct x of the mesh points (wedge along x;
    the cell centroids are not the midpoints)."""
    c = values(os.path.join(case, time, "Ccx"), 0)[1]
    s = open(os.path.join(case, "constant", "polyMesh", "points")).read()
    pts = np.array(re.findall(r"\(([^()]*)\)", s[s.index("(", s.index("}")):]),
                   dtype=str)
    x = np.array([float(p.split()[0]) for p in pts])
    e = np.unique(np.round(x, 12))
    assert len(e) == len(c) + 1, (len(e), len(c))
    assert np.all((c > e[:-1]) & (c < e[1:]))
    return c, e


cs, es = edges(src, st)
ct, et = edges(tgt, tt)
ns, nt = len(cs), len(ct)
print(f"source {ns} cells to {es[-1]:.4f} m, target {nt} cells to {et[-1]:.4f} m")

# overlap matrix W[t, s] in volume (r^3) units
vs = es[1:]**3 - es[:-1]**3
W = np.zeros((nt, ns))
for i in range(nt):
    lo, hi = et[i], et[i + 1]
    j0 = max(np.searchsorted(es, lo, "right") - 1, 0)
    j1 = min(np.searchsorted(es, hi, "left"), ns)
    for j in range(j0, j1):
        a, b = max(lo, es[j]), min(hi, es[j + 1])
        if b > a:
            W[i, j] = b**3 - a**3
    if hi > es[-1]:                       # beyond the source: outermost cell
        W[i, ns - 1] += hi**3 - max(lo, es[-1])**3
vt = W.sum(axis=1)
assert np.allclose(vt, et[1:]**3 - et[:-1]**3, rtol=1e-9)

sdir = os.path.join(src, st)
fields = [f for f in os.listdir(sdir)
          if os.path.isfile(os.path.join(sdir, f)) and f not in ("Ccx", "Ccy", "Ccz", "C", "V", "Vc")
          and not f.startswith("alphaRho") and not f.startswith("phi")]
S = {f: values(os.path.join(sdir, f), ns) for f in fields}
phases = sorted({f.split(".", 1)[1] for f in fields if f.startswith("alpha.")})


def vavg(v, w=None):
    w = np.ones(ns) if w is None else w
    num = W @ (v*w[:, None] if v.ndim == 2 else v*w)
    den = W @ w
    out_v = W @ v / vt[:, None] if v.ndim == 2 else W @ v / vt
    ok = den > 1e-300*vt
    if v.ndim == 2:
        res = out_v.copy()
        res[ok] = num[ok]/den[ok, None]
    else:
        res = out_v.copy()
        res[ok] = num[ok]/den[ok]
    return res


out = {}
for f in fields:
    kind, v = S[f]
    phase = f.split(".", 1)[1] if "." in f else None
    if f == "pTemporalMax":
        out[f] = np.array([v[W[i] > 0].max() for i in range(nt)])
    elif f.startswith("alpha.") or phase is None and f != "U":
        out[f] = vavg(v)
    elif f == "U":
        out[f] = vavg(v, S["rho"][1])
    elif phase in phases and f == f"rho.{phase}":
        out[f] = vavg(v, S[f"alpha.{phase}"][1])
    elif phase in phases:
        out[f] = vavg(v, S[f"alpha.{phase}"][1]*S[f"rho.{phase}"][1])
    else:
        out[f] = vavg(v)
    tdir = os.path.join(tgt, tt)
    write(os.path.join(tdir, f), os.path.join(sdir, f), kind, out[f])

# checks: total mass and momentum of each phase
for ph in phases:
    a, r = S[f"alpha.{ph}"][1], S[f"rho.{ph}"][1]
    m_s = (W.sum(axis=0)*a*r).sum()      # source incl. the extension
    m_t = (vt*out[f"alpha.{ph}"]*out[f"rho.{ph}"]).sum()
    print(f"{ph}: mass source {m_s:.6e} target {m_t:.6e}")
sa = sum(out[f"alpha.{ph}"] for ph in phases)
print(f"sum alpha: {sa.min():.12f} - {sa.max():.12f}; mapped {len(fields)} fields")
