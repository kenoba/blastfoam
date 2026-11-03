#!/usr/bin/env python3
"""
compare_pontalier.py: arrival times at 1.0, 1.5, 2.0 m and the 2.0 m peak
overpressure of Pontalier cases against the measurements (pontalier_data.py).

    python3 compare_pontalier.py <case> [<case> ...]

Arrival: first 5% overpressure at the probe (the paper takes the time of the
maximum overpressure; on the graded mesh the model peak lags the front by a
few cells, so both are printed). Ratios are relative to the baseline case
(trial 0) when it is among the cases, as in the paper's normalised plots.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pontalier_data import DATA, GAUGES  # noqa: E402


def read_case(case):
    info = open(os.path.join(case, "charge.txt")).read().split()
    kv = dict(zip(info[0::2], info[1::2]))
    p0 = float(kv["p0"])
    path = os.path.join(case, "postProcessing", "pressureProbes", "0", "p")
    radii = [float(l.split("(")[1].split()[0]) for l in open(path) if l.startswith("# Probe")]
    d = np.loadtxt(path)
    t = d[:, 0]
    res = {}
    for i, r in enumerate(radii):
        dp = d[:, i + 1] - p0
        hit = dp > 0.05*p0
        if not hit.any():
            continue
        k = np.argmax(hit)
        # first local maximum after arrival: the shock peak
        j = k + int(np.argmax(dp[k:k + max(1, int(np.searchsorted(t, t[k] + 2e-4) - k))]))
        ended = (dp[k:] < 0).any()
        neg = np.argmax(dp[k:] < 0) if ended else len(dp) - k
        imp = np.trapezoid(np.clip(dp[k:k + neg], 0, None), t[k:k + neg])
        res[round(r, 3)] = dict(ta=t[k], tpk=t[j], pk=dp.max(), imp=imp if ended else np.nan)
    return kv, res


def main():
    cases = sys.argv[1:]
    if not cases:
        sys.exit(__doc__)
    rows = []
    for c in cases:
        kv, res = read_case(c)
        rows.append((os.path.basename(os.path.normpath(c)), kv, res))
    print(f"{'case':<12}{'trial':>6}{'closure':>9}  "
          + "".join(f"{'ta ' + str(g) + ' m':>16}" for g in GAUGES)
          + f"{'peak 2 m [kPa]':>26}{'I+ 2 m':>9}")
    print(f"{'':<27}" + "".join(f"{'model/meas [ms]':>16}" for _ in GAUGES)
          + f"{'model/gauge/video':>26}{'[Pa s]':>9}")
    for name, kv, res in rows:
        tr = kv["trial"]
        meas = DATA[tr]
        line = f"{name:<12}{tr:>6}{kv['closure']:>9}  "
        for g, tm in zip(GAUGES, meas[0]):
            r = res.get(g)
            line += (f"{r['ta']*1e3:7.3f}/{tm if tm else float('nan'):6.3f}  " if r else f"{'-':>16}")
        r = res.get(2.0)
        gp = meas[1] if meas[1] else float("nan")
        line += (f"{r['pk']/1e3:10.1f}/{gp:6.1f}/{meas[2]:6.1f}  {r['imp']:7.1f}" if r else "")
        print(line)
    # normalised by the baseline of the same closure
    base = {kv["closure"]: res for _, kv, res in rows if kv["trial"] == "0"}
    if base:
        print("\nnormalised by trial 0 (model | measured): ta at 2.0 m, peak at 2.0 m (gauge)")
        b0 = DATA["0"]
        for name, kv, res in rows:
            b = base.get(kv["closure"])
            if not b or kv["trial"] == "0" or 2.0 not in res or 2.0 not in b:
                continue
            m = DATA[kv["trial"]]
            ta_m = res[2.0]["ta"]/b[2.0]["ta"]
            pk_m = res[2.0]["pk"]/b[2.0]["pk"]
            ta_e = m[0][2]/b0[0][2] if m[0][2] else float("nan")
            pk_e = m[1]/b0[1] if m[1] else float("nan")
            print(f"{name:<12} ta {ta_m:6.3f} | {ta_e:6.3f}   peak {pk_m:6.3f} | {pk_e:6.3f}")


if __name__ == "__main__":
    main()
