#!/usr/bin/env python3
"""First-positive-phase metrics of the 1D mixing runs (mix/c*) against the
laminar baseline s4_lam_new: arrival, peak, duration and impulse of the first
positive overpressure phase at every probe radius.
Balakrishnan et al. (2010) give the impulse of the first positive phase;
monitor_mix.py integrates all positive overpressure, later pulses included.

usage: FOI2013_RUNS=<run dir> first_phase.py [case ...]   (default mix/c*)"""
import glob
import os
import re
import sys

import numpy as np
if not hasattr(np, "trapezoid"):  # numpy < 2.0
    np.trapezoid = np.trapz

ROOT = os.environ.get("FOI2013_RUNS", os.getcwd())
p0 = 101325.0


def probes(case):
    base = os.path.join(case, "postProcessing", "pressureProbes")
    parts, radii = [], None
    for start in sorted(os.listdir(base), key=float):
        path = os.path.join(base, start, "p")
        radii = [float(re.search(r"\(([^ ]+)", l).group(1)) for l in open(path) if l.startswith("# Probe")]
        a = np.loadtxt(path, ndmin=2, comments="#")
        parts.append(a[a[:, 0] > parts[-1][-1, 0]] if parts else a)
    return radii, np.vstack(parts)


def secondary_shock(t, dp, i0, rise=10e3, window=10e-6, delay=50e-6):
    """Index of the secondary shock: first rise of more than `rise` within
    `window`, later than `delay` after the main peak (the peak within 0.2 ms
    of the arrival at i0; None if none)."""
    ip = i0 + int(np.argmax(dp[i0:np.searchsorted(t, t[i0] + 2e-4)]))
    for j in range(ip, len(t)):
        if t[j] < t[ip] + delay:
            continue
        k = np.searchsorted(t, t[j] - window)
        if dp[j] - dp[k:j + 1].min() > rise:
            return k + int(np.argmin(dp[k:j + 1]))
    return None


def first_phase(t, dp, thr=0.05):
    """Arrival: first dp > thr*p0. Phase: from the last dp <= 0 before it to
    the first dp < 0 after it or the foot of the secondary shock, whichever
    comes first (the mixing runs stay positive until the secondary shock at
    0.75 m); open if the record ends first."""
    on = np.nonzero(dp > thr * p0)[0]
    if len(on) == 0:
        return None
    i0 = on[0]
    neg = np.nonzero(dp[:i0] <= 0)[0]
    s = neg[-1] + 1 if len(neg) else 0
    off = np.nonzero(dp[i0:] < 0)[0]
    e = i0 + off[0] if len(off) else len(dp)
    i2 = secondary_shock(t, dp, i0)
    cut = i2 is not None and i2 < e
    if cut:
        e = i2 + 1
    return dict(ta=t[i0], pk=dp[s:e].max(), tend=t[e - 1], open=not len(off) and not cut,
                cut=cut, I=np.trapezoid(dp[s:e], t[s:e]))


def table(case):
    radii, d = probes(case)
    return {round(r, 3): first_phase(d[:, 0], d[:, i + 1] - p0) for i, r in enumerate(radii)}


if __name__ == "__main__":
    cases = sys.argv[1:] or sorted(glob.glob(os.path.join(ROOT, "mix", "c*")))
    base = table(os.path.join(ROOT, "s4_lam_new"))
    for case in cases:
        print(os.path.basename(case.rstrip("/")))
        for r, m in table(case).items():
            if m is None:
                continue
            b = base.get(r)
            rel = (f"  peak {(m['pk'] / b['pk'] - 1) * 100:+6.1f}%  I {(m['I'] / b['I'] - 1) * 100:+6.1f}%"
                   f"{' (baseline open)' if b['open'] else ''}") if b else "  (no baseline probe)"
            print(f"   r {r:5.3f} m: arrival {m['ta'] * 1e3:6.3f} ms, end {m['tend'] * 1e3:6.3f} ms"
                  f"{' OPEN' if m['open'] else (' 2ND ' if m['cut'] else '     ')}, peak {m['pk'] / 1e3:7.1f} kPa, I {m['I']:6.1f} Pa s{rel}")
    print("s4_lam_new")
    for r, m in base.items():
        if m:
            print(f"   r {r:5.3f} m: arrival {m['ta'] * 1e3:6.3f} ms, end {m['tend'] * 1e3:6.3f} ms"
                  f"{' OPEN' if m['open'] else (' 2ND ' if m['cut'] else '     ')}, peak {m['pk'] / 1e3:7.1f} kPa, I {m['I']:6.1f} Pa s")
