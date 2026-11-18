#!/usr/bin/env python3
"""Ground pressure at P2/P3 (4/6 m): FOI shot 4 measurement against the 2D run.
    waveforms.py <2D case> <out.png>"""
import sys, glob, re, os
import numpy as np
if not hasattr(np, "trapezoid"):  # numpy < 2.0
    np.trapezoid = np.trapz
sys.path.insert(0, os.environ.get("BLAST_DIR", "/opt/OpenFOAM/blastfoam") + "/tutorials/blastEulerFoam/FOI2013_TNT_1D_spherical")
import extract_shots as ex
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
case, out = sys.argv[1:3]
p0 = 101325.0
tm, pm, _ = ex.load(4)
fs = sorted(glob.glob(f"{case}/postProcessing/pressureProbes/*/p"), key=lambda f: float(f.split('/')[-2]))
xs = [float(re.search(r'\(([-\d.e]+)', l).group(1)) for l in open(fs[0]) if l.startswith('# Probe')]
d = np.vstack([np.loadtxt(f, comments='#', ndmin=2) for f in fs]); d = d[np.argsort(d[:, 0], kind='stable')]
_, u = np.unique(d[:, 0], return_index=True); d = d[u]
def feats(t, p, thr):
    ka = np.where((t > 0.2e-3) & (p > thr))[0][0]
    while ka > 0 and p[ka-1] > 0.25*thr: ka -= 1
    ta = t[ka]
    dt = np.median(np.diff(t)); w = max(int(0.2e-3/dt), 1)
    ps = np.convolve(p, np.ones(w)/w, "same")
    j = ka + max(int(0.05e-3/dt), 1); j += int(np.argmax(ps[j:] < 0))
    I = np.trapezoid(p[ka:j], t[ka:j])
    pk = p[(t >= ta) & (t < ta + 0.3e-3)].max()
    # cumulative impulse at fixed times after arrival
    cum = [np.trapezoid(p[ka:][(t[ka:]-ta) <= s], t[ka:][(t[ka:]-ta) <= s]) for s in (0.25e-3, 0.5e-3, 1e-3, 2e-3)]
    # negative phase: min and time of the secondary shock (largest positive jump after t+)
    neg = (t > t[j]) & (t < t[j] + 8e-3)
    pmin = ps[neg].min(); tmin = t[neg][np.argmin(ps[neg])] - ta
    return dict(ta=ta, pk=pk, td=t[j]-ta, I=I, cum=cum, pmin=pmin, tmin=tmin, ka=ka)
fig, ax = plt.subplots(2, 1, figsize=(9, 7))
for k, (r, thr) in enumerate(((4.0, 8e3), (6.0, 4e3))):
    tmod, pmod = d[:, 0], d[:, 1 + xs.index(r)] - p0
    a = feats(tm, pm[:, k+1], thr); b = feats(tmod, pmod, thr)
    print(f"P{k+2} ({r:.0f} m)       ta [ms]  peak [kPa]  t+ [ms]  I+ [Pa s]  I(0.25/0.5/1/2 ms) [Pa s]          neg min [kPa] at [ms]")
    for n, f in (("measured", a), ("2D model", b)):
        print(f"  {n:9s} {f['ta']*1e3:8.3f} {f['pk']/1e3:10.1f} {f['td']*1e3:8.2f} {f['I']:9.1f}  "
              + " ".join(f"{c:6.1f}" for c in f['cum']) + f"   {f['pmin']/1e3:8.1f} {f['tmin']*1e3:6.2f}")
    print("  ratio    " + f"{'':8s} {b['pk']/a['pk']:10.2f} {b['td']/a['td']:8.2f} {b['I']/a['I']:9.2f}  "
          + " ".join(f"{y/x:6.2f}" for x, y in zip(a['cum'], b['cum'])))
    ax[k].plot((tm - a['ta'])*1e3, pm[:, k+1]/1e3, lw=0.6, label="measured")
    ax[k].plot((tmod - b['ta'])*1e3, pmod/1e3, lw=1.2, label="2D model")
    ax[k].set_xlim(-0.5, 12); ax[k].axhline(0, c='k', lw=0.4); ax[k].set_ylabel(f"P{k+2} overpressure [kPa]")
    ax[k].legend()
ax[1].set_xlabel("time after arrival [ms]")
fig.tight_layout(); fig.savefig(out, dpi=120)
