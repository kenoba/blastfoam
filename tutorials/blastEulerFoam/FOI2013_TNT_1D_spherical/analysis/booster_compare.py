#!/usr/bin/env python3
"""Booster (top) vs centre (ctr) initiation in the 2D free field: positive-phase
impulse, hump times after arrival, and traces at the gauge angles."""
import os, sys, importlib.util, numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
spec = importlib.util.spec_from_file_location("m", os.path.join(os.path.dirname(os.path.abspath(__file__)), "monitor_free2d.py"))
src = open(spec.origin).read().split("res = {}")[0]; m = {}; exec(src, m)
P0 = 101325.0
def key(p):
    r = np.hypot(*p); a = np.degrees(np.arctan2(p[1], p[0])); return round(r, 3), round(a, 1)
data = {}
for c in ("ctr", "top"):
    pos, d = m["probes"](c); data[c] = {key(p): (d[:, 0], d[:, i+1]) for i, p in enumerate(pos)}
def analyse(t, p):
    dp = p - P0; ka = int(np.argmax(dp > 0.05*dp.max())); ta = t[ka]
    kn = ka + int(np.argmax(dp[ka:] < 0)) if (dp[ka:] < 0).any() else len(t)-1
    Ip = np.trapezoid(dp[ka:kn], t[ka:kn])
    msk = (t > ta+0.2e-3) & (t < ta+2.5e-3); tt, pp = t[msk], dp[msk]
    g = np.linspace(tt[0], tt[-1], 500); s = np.convolve(np.interp(g, tt, pp), np.ones(10)/10, "same")
    hs = [(g[i]-ta)*1e3 for i in range(1, len(s)-1) if s[i] > s[i-1] and s[i] >= s[i+1] and s[i]-s[max(0,i-50):i].min() > 0.03*dp.max()]
    return ta, Ip, (t[kn]-ta)*1e3, hs
print(f"{'r':>6} {'ang':>6} | {'I+ ctr':>7} {'I+ top':>7} {'dI+':>6} | {'t+ ctr':>6} {'t+ top':>6} | humps ctr / top (ms after arrival)")
for (r, a) in sorted(data["ctr"]):
    if r < 1.0: continue
    A = analyse(*data["ctr"][(r, a)]); B = analyse(*data["top"][(r, a)])
    if data["ctr"][(r,a)][0][-1] < A[0] + 2.5e-3: hs = ""
    else: hs = f"{[round(h,2) for h in A[3]]} / {[round(h,2) for h in B[3]]}"
    print(f"{r:6.3f} {a:6.1f} | {A[1]:7.1f} {B[1]:7.1f} {100*(B[1]/A[1]-1):+5.1f}% | {A[2]:6.2f} {B[2]:6.2f} | {hs}")
fig, ax = plt.subplots(2, 3, figsize=(14, 7), sharex=True)
for j, r in enumerate((1.0, 1.5, 2.0)):
    for i, a in enumerate((-18.4, -45.0)):
        for c, ls in (("ctr", "-"), ("top", "--")):
            t, p = data[c][(r, a)]; ta = analyse(t, p)[0]
            ax[i, j].plot((t-ta)*1e3, (p-P0)/1e3, ls, label=c)
        ax[i, j].set_title(f"r {r} m, {a} deg"); ax[i, j].set_xlim(-0.1, 3); ax[i, j].grid(alpha=.3)
        for h in (0.6, 0.9, 2.0, 2.3): ax[i, j].axvline(h, color="k", lw=.5, alpha=.4)
ax[0, 0].legend(); [x.set_xlabel("t - t_arrival [ms]") for x in ax[1]]; ax[0,0].set_ylabel("dp [kPa]"); ax[1,0].set_ylabel("dp [kPa]")
fig.suptitle("2D free field: centre vs booster initiation (thin lines: measured hump times)")
fig.tight_layout(); fig.savefig(sys.argv[1], dpi=110)
