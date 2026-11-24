#!/usr/bin/env python3
"""Ground overpressure at P1-P3 (2/4/6 m), FOI shots 4 (TNT) and 5 (TNT/Al):
measurement against two 2D runs.
    waveforms_al.py <2D shot 4 case> <2D shot 5 case> <out.pdf>"""
import sys, glob, re, os
import numpy as np
sys.path.insert(0, os.environ.get("BLAST_DIR", "/opt/OpenFOAM/blastfoam")
                + "/tutorials/blastEulerFoam/FOI2013_TNT_1D_spherical")
import extract_shots as ex
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

c4, c5, out = sys.argv[1:4]
p0 = 101325.0


def model(case):
    fs = sorted(glob.glob(f"{case}/postProcessing/pressureProbes/*/p"),
                key=lambda f: float(f.split('/')[-2]))
    xs = [(float(m.group(1)), float(m.group(2))) for l in open(fs[0]) if l.startswith('# Probe')
          for m in [re.search(r'\(([-\d.e]+) ([-\d.e]+)', l)]]
    d = np.vstack([np.loadtxt(f, comments='#', ndmin=2) for f in fs])
    d = d[np.argsort(d[:, 0], kind='stable')]
    _, u = np.unique(d[:, 0], return_index=True)
    d = d[u]
    return d[:, 0], {r: d[:, 1 + i] - p0 for i, (r, y) in enumerate(xs) if y < 0.01}


def smooth(p, n):
    return np.convolve(p, np.ones(n)/n, "same")


fig, ax = plt.subplots(3, 1, figsize=(7.5, 9))
t4, m4 = model(c4)
t5, m5 = model(c5)
for k, (r, win, ylim) in enumerate(((2.0, (3.0, 9.0), (-30, 350)), (4.0, (6.8, 13.5), (-15, 130)),
                                    (6.0, (11.5, 18.0), (-10, 55)))):
    a = ax[k]
    for shot, col, lab in ((4, "0.6", "shot 4 TNT, measured"), (5, "k", "shot 5 TNT/Al, measured")):
        tm, pm, _ = ex.load(shot)
        a.plot(tm*1e3, smooth(pm[:, k], 20)/1e3, color=col, lw=0.8, label=lab)
    a.plot(t4*1e3, m4[r]/1e3, color="C0", lw=1.2, ls="--", label="shot 4, model")
    a.plot(t5*1e3, m5[r]/1e3, color="C3", lw=1.2, label="shot 5, model")
    a.set_xlim(*win)
    a.set_ylim(*ylim)
    a.axhline(0, color="0.8", lw=0.5)
    a.set_ylabel("overpressure [kPa]")
    a.set_title(f"P{k+1} ({r:.0f} m)", fontsize=10)
ax[-1].set_xlabel("time after detonation [ms]")
ax[0].legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig(out)
print("wrote", out)
