#!/usr/bin/env python3
"""Seed a perturbation of the products/air contact of the mapped 2D state.

The 2D fields at t0 are the 1D spherical solution mapped by rotateFields.
In a band around the contact radius r_c (|r - r_c| < h) every cell is
reset to the 1D state at the displaced radius

    r' = r - delta(psi) w(r),   w(r) = cos^2(pi (r - r_c)/(2 h)),

where psi is the angle from the symmetry plane (y = 0) and delta(psi) a sum
of cos(n psi) modes with even n (symmetric about the plane and smooth on
the axis). The state is taken from the 1D cell containing r' (piecewise
constant), so every cell holds a consistent 1D state (primitive and
conserved fields together); vectors are the 1D radial component along the
2D radial direction. The lead shock and the secondary shock must lie
outside the band: the script stops if the pressure jump of either is
within it.

With "charge" the perturbation is instead put on the surface of the charge
at t = 0 (run from 0): setRefinedFields sets the charge as a sphere of
radius r0 + max(delta) at the finest level, and the cells with
r > r0 + delta(psi) are reset to the air state (taken from the outermost
cell); the surface is sharp, as set by setFields.

Usage (in the case directory, on the decomposed case with cell centres C
written by postProcess -func writeCellCentres):
  perturb_contact.py map <1D time dir> <t0> <seed> [--halfwidth h]
  perturb_contact.py charge <r0> <seed>
  perturb_contact.py maxdisp <seed>        (prints max |delta| [m])
<seed>: --modes n:amp [n:amp ...] (amplitudes in m) or
        --multimode nmin nmax rms [--seed s] (all even n in [nmin, nmax]
        with normal random amplitudes scaled to the rms displacement).
The contact radius is where the N2 mass fraction of the gas is halfway
between its products and air values; h defaults to 8 mm.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

HDR = re.compile(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*(\d+)\s*\(")
UNI = re.compile(r"internalField\s+uniform\s+([^;]+);")


def read_field(path):
    """Return (kind, values as list of str tokens or None, text, span)."""
    text = open(path).read()
    m = HDR.search(text)
    if m:
        kind, n = m.group(1), int(m.group(2))
        end = text.index("\n)", m.end())
        body = text[m.end():end]
        if kind == "scalar":
            vals = body.split()
        else:
            vals = re.findall(r"\(([^()]*)\)", body)
        if len(vals) != n:
            sys.exit(f"{path}: expected {n} values, read {len(vals)}")
        return kind, vals, text, (m.start(), end + 2)
    u = UNI.search(text)
    if u:
        v = u.group(1).strip()
        kind = "vector" if v.startswith("(") else "scalar"
        return kind, v, text, (u.start(), u.end())
    return None, None, text, None


def write_field(path, kind, vals, text, span):
    if kind == "scalar":
        body = "\n".join(vals)
    else:
        body = "\n".join(f"({v})" for v in vals)
    new = f"internalField   nonuniform List<{kind}> \n{len(vals)}\n(\n{body}\n)"
    # A uniform entry's span includes its ';'
    if text[span[1] - 1] == ";":
        new += ";"
    with open(path, "w") as f:
        f.write(text[:span[0]] + new + text[span[1]:])


def centres(path):
    kind, vals, _, _ = read_field(path)
    return np.array([list(map(float, v.split())) for v in vals])


def parse_modes(a):
    modes = []
    for m in a.modes:
        n, amp = m.split(":")
        modes.append((int(n), float(amp)))
    if a.multimode:
        nmin, nmax, rms = int(a.multimode[0]), int(a.multimode[1]), a.multimode[2]
        ns = [n for n in range(nmin, nmax + 1) if n % 2 == 0]
        rng = np.random.default_rng(a.seed)
        amps = rng.normal(size=len(ns))
        amps *= rms/np.sqrt(0.5*np.sum(amps**2))
        modes += list(zip(ns, amps))
    if not modes or any(n % 2 for n, _ in modes):
        sys.exit("need modes with even n")
    return modes


def cells(tdir):
    C = centres(os.path.join(tdir, "C"))
    r = np.linalg.norm(C, axis=1)
    psi = np.arctan2(C[:, 1], np.hypot(C[:, 0], C[:, 2]))
    return C, r, psi


def seed_charge(a, delta, dmax):
    """Reset cells beyond r0 + delta(psi) of the enlarged charge to air."""
    r0 = a.r0
    total = 0
    procs = sorted(glob.glob("processor*")) or ["."]
    # Air state: the outermost cell of the whole mesh
    rmax, air = -1.0, None
    for pd in procs:
        C, r, psi = cells(os.path.join(pd, "0"))
        k = int(np.argmax(r))
        if r[k] > rmax:
            rmax, air = r[k], (pd, k)
    for pd in procs:
        tdir = os.path.join(pd, "0")
        C, r, psi = cells(tdir)
        sel = np.where((r < r0 + dmax + 1e-3) & (r > r0 + delta(psi)))[0]
        total += len(sel)
        if not len(sel):
            continue
        for f in os.listdir(tdir):
            path = os.path.join(tdir, f)
            if not os.path.isfile(path) or f in ("C", "Ccx", "Ccy", "Ccz"):
                continue
            kind, vals, text, span = read_field(path)
            if kind is None or span is None:
                continue
            akind, avals, _, _ = read_field(os.path.join(air[0], "0", f))
            aval = avals[air[1]] if isinstance(avals, list) else \
                (avals.strip("()") if akind == "vector" else avals)
            if not isinstance(vals, list):
                if vals == (f"({aval})" if kind == "vector" else aval):
                    continue
                vals = [vals.strip("()") if kind == "vector" else vals]*len(r)
            for c in sel:
                vals[c] = aval
            write_field(path, kind, vals, text, span)
    print(f"charge radius {r0} m + perturbation (max {dmax*1e3:.2f} mm): "
          f"reset {total} cells to air in {len(procs)} processor directories")


def seed_map(a, delta, dmax):
    """Displace the contact of the mapped 1D state."""
    x1 = centres(os.path.join(a.oneD, "C"))[:, 0]
    faces = np.concatenate(([0.0], 0.5*(x1[1:] + x1[:-1]), [np.inf]))
    n2 = np.array(read_field(os.path.join(a.oneD, "N2.gas"))[1], float)
    p1 = np.array(read_field(os.path.join(a.oneD, "p"))[1], float)
    # Products: the N2 fraction where the gas holds most CO; air: the outer
    # boundary. Search inwards from the outside (inside the charge the gas
    # is a trace with its own composition).
    co = np.array(read_field(os.path.join(a.oneD, "CO.gas"))[1], float)
    n2p, n2a = np.median(n2[co > 0.9*co.max()]), n2[-1]
    half = 0.5*(n2p + n2a)
    k = len(n2) - 1
    while n2[k - 1] > half:
        k -= 1
    k -= 1     # n2[k] <= half < n2[k + 1]
    rc = x1[k] + (half - n2[k])/(n2[k + 1] - n2[k])*(x1[k + 1] - x1[k])
    h = a.halfwidth
    if dmax*np.pi/(2*h) >= 0.9:
        sys.exit(f"displacement {dmax} too large for the band half-width {h}")
    # Pressure jumps (> 30% between neighbouring cells) must lie outside the
    # band (the displaced radius r' stays within it: w = 0 at its edges)
    band = (x1 > rc - h) & (x1 < rc + h)
    jump = np.abs(np.diff(p1))/np.minimum(p1[1:], p1[:-1]) > 0.3
    bad = jump & band[1:] & band[:-1]
    if np.any(bad):
        sys.exit(f"pressure jump at r = {x1[1:][bad]} m within the band "
                 f"r_c +- h = {rc:.4f} +- {h} m")
    print(f"contact radius {rc:.5f} m, band +-{h} m, max displacement "
          f"{dmax*1e3:.2f} mm")

    fields1 = {}
    for f in os.listdir(a.oneD):
        p = os.path.join(a.oneD, f)
        if os.path.isfile(p):
            kind, vals, _, _ = read_field(p)
            if kind and isinstance(vals, list) and len(vals) == len(x1):
                fields1[f] = (kind, vals)

    procs = sorted(glob.glob("processor*")) or ["."]
    total = 0
    for pd in procs:
        tdir = os.path.join(pd, a.t0)
        C, r, psi = cells(tdir)
        sel = np.where(np.abs(r - rc) < h)[0]
        total += len(sel)
        if not len(sel):
            continue
        rs = r[sel]
        w = np.cos(0.5*np.pi*(rs - rc)/h)**2
        r1 = rs - delta(psi[sel])*w
        idx = np.searchsorted(faces, r1) - 1
        er = C[sel]/rs[:, None]
        for f, (kind1, vals1) in fields1.items():
            path = os.path.join(tdir, f)
            if not os.path.exists(path) or f in ("C", "Ccx", "Ccy", "Ccz"):
                continue
            kind, vals, text, span = read_field(path)
            if kind != kind1 or span is None:
                continue
            if not isinstance(vals, list):
                vals = [vals.strip("()") if kind == "vector" else vals]*len(r)
            if kind == "scalar":
                for c, i in zip(sel, idx):
                    vals[c] = vals1[i]
            else:
                ur = np.array([float(vals1[i].split()[0]) for i in idx])
                for j, c in enumerate(sel):
                    u = ur[j]*er[j]
                    vals[c] = f"{u[0]:.6g} {u[1]:.6g} {u[2]:.6g}"
            write_field(path, kind, vals, text, span)
    print(f"reset {total} cells in {len(procs)} processor directories")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    pm = sub.add_parser("map")
    pm.add_argument("oneD")
    pm.add_argument("t0")
    pm.add_argument("--halfwidth", type=float, default=0.008)
    pc = sub.add_parser("charge")
    pc.add_argument("r0", type=float)
    pd = sub.add_parser("maxdisp")
    for p in (pm, pc, pd):
        p.add_argument("--modes", nargs="*", default=[])
        p.add_argument("--multimode", nargs=3, type=float)
        p.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    modes = parse_modes(a)

    def delta(psi):
        return sum(amp*np.cos(n*psi) for n, amp in modes)

    dmax = max(abs(delta(np.linspace(0, np.pi/2, 2001))))
    if a.cmd == "maxdisp":
        print(f"{dmax:.6g}")
    elif a.cmd == "charge":
        seed_charge(a, delta, dmax)
    else:
        seed_map(a, delta, dmax)


if __name__ == "__main__":
    main()
