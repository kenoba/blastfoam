#!/usr/bin/env python3
"""Monitor of the 1D mixing runs (mix/c*): health, burnt fraction from O2
(o2burn.py at the latest write) and probe peak/impulse against the laminar
baseline s4_lam_new at 0.75, 1.0 and 1.5 m."""
import glob
import os
import re
import subprocess
import time

import numpy as np
if not hasattr(np, "trapezoid"):  # numpy < 2.0
    np.trapezoid = np.trapz

ROOT = os.environ.get("FOI2013_RUNS", os.getcwd())  # directory with mix/ and s4_lam_new/
ENV = ("source $WM_PROJECT_DIR/etc/bashrc >/dev/null 2>&1; "
       "source $BLAST_DIR/etc/bashrc >/dev/null 2>&1; ")
p0 = 101325.0


def probes(case):
    base = os.path.join(case, "postProcessing", "pressureProbes")
    parts = []
    radii = None
    for start in sorted(os.listdir(base), key=float):
        path = os.path.join(base, start, "p")
        radii = [float(re.search(r"\(([^ ]+)", l).group(1)) for l in open(path) if l.startswith("# Probe")]
        a = np.loadtxt(path, ndmin=2, comments="#")
        parts.append(a[a[:, 0] > parts[-1][-1, 0]] if parts else a)
    return radii, np.vstack(parts)


def metrics(case, r, tmax):
    radii, d = probes(case)
    i = int(np.argmin([abs(x - r) for x in radii]))
    t = d[:, 0]
    m = t <= tmax
    dp = d[m, i + 1] - p0
    if dp.max() < 0.05 * p0:
        return None
    pk = dp.max() / 1e3
    I = np.trapezoid(np.clip(dp, 0, None), t[m])
    return pk, I


def latest(case):
    ts = [d for d in os.listdir(case) if re.fullmatch(r"[0-9.e-]+", d) and d != "0"]
    return max(ts, key=float) if ts else None


base = os.path.join(ROOT, "s4_lam_new")
for case in sorted(glob.glob(os.path.join(ROOT, "mix", "c*"))):
    name = os.path.basename(case)
    log = os.path.join(case, "log.blastEulerFoam")
    if not os.path.isfile(log):
        print(f"{name}: no log")
        continue
    txt = open(log, errors="ignore").read()[-200000:]
    ts = re.findall(r"\nTime = ([0-9.e-]+)", txt)
    dts = re.findall(r"deltaT = ([0-9.e-]+)", txt)
    et = re.findall(r"ExecutionTime = ([0-9.]+) s", txt)
    tnow = float(ts[-1]) if ts else 0.0
    age = time.time() - os.path.getmtime(log)
    bad = [k for k in ("FOAM FATAL", "Floating point", "nan", "NaN") if k in txt]
    done = "End" in txt[-500:]
    print(f"{name}: t = {tnow * 1e3:.3f} ms, dt = {dts[-1] if dts else '?'}, CPU {float(et[-1]) / 3600 if et else 0:.2f} h, "
          f"log age {age:.0f} s{', ERRORS ' + str(bad) if bad else ''}{', finished' if done else ''}")
    lt = latest(case)
    if lt:
        subprocess.run(["bash", "-c", ENV + f"cd {case} && postProcess -func writeCellVolumes -time {lt} > log.vol 2>&1"])
        if os.path.isdir(os.path.join(base, lt)) and not os.path.isfile(os.path.join(base, lt, "Vc")):
            subprocess.run(["bash", "-c", ENV + f"cd {base} && postProcess -func writeCellVolumes -time {lt} > log.vol 2>&1"])
        out = subprocess.run(["python3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "o2burn.py"), case, lt],
                             capture_output=True, text=True).stdout.strip().splitlines()
        bout = subprocess.run(["python3", os.path.join(os.path.dirname(os.path.abspath(__file__)), "o2burn.py"), base, lt],
                              capture_output=True, text=True).stdout.strip().splitlines()
        print(f"   burnt at {float(lt) * 1e3:.3f} ms: {out[-1].split()[-1] if len(out) > 1 else '?'}"
              f" (baseline {bout[-1].split()[-1] if len(bout) > 1 else 'n/a'})")
    for r in (0.75, 1.0, 1.5):
        a = metrics(case, r, tnow)
        b = metrics(base, r, tnow)
        if a and b:
            print(f"   r {r:4.2f} m: peak {a[0]:7.1f} kPa ({(a[0] / b[0] - 1) * 100:+5.1f}%), "
                  f"I(t<={tnow * 1e3:.2f} ms) {a[1]:6.1f} Pa s ({(a[1] / b[1] - 1) * 100:+5.1f}%)")
