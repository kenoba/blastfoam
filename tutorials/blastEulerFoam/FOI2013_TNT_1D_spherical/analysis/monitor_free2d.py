#!/usr/bin/env python3
"""Status and physics check of the 2D free-field booster runs (free2d/ctr, top).
Run health (errors, process, dt, rate, ETA, cells) and the case physics:
ctr against the 1D run (arrival) and isotropy over the probe angles; top
against ctr per angle (arrival, peak, impulse in the first 1.5 ms, humps).
Prints flags that call for stopping or changing the runs."""
import os, re, subprocess, time, json, glob
import numpy as np
D = os.path.join(os.environ.get("FOI2013_RUNS", os.getcwd()), "free2d")
STATE = os.path.join(D, "monitor_state.json")
P0 = 101325.0
TA1D = {0.5: 0.154e-3, 1.0: 0.515e-3, 1.5: 1.122e-3, 2.0: 1.951e-3, 2.5: 2.933e-3, 2.828: 3.619e-3}
st = json.load(open(STATE)) if os.path.exists(STATE) else {}
now = time.time(); flags = []; out = []

def tail_val(log, pat):
    r = subprocess.run(["bash", "-c", f"grep -a -E '{pat}' {log} | tail -1"], capture_output=True, text=True).stdout
    return r.strip()

def probes(c):
    fs = sorted(glob.glob(os.path.join(D, c, "postProcessing/pressureProbes/*/p")),
                key=lambda f: float(f.split("/")[-2]))
    lab = [l for l in open(fs[0]) if l.startswith("# Probe")]
    pos = [tuple(map(float, re.search(r"\(([^)]*)\)", l).group(1).split()[:2])) for l in lab]
    parts = [np.loadtxt(f, comments="#", ndmin=2) for f in fs]
    # a restart overrides the earlier file from its start time on
    d = parts[0]
    for f, a in zip(fs[1:], parts[1:]):
        d = np.vstack([d[d[:, 0] < float(f.split("/")[-2])], a])
    return pos, d

def feats(t, p):
    dp = p - P0
    if dp.max() < 0.05*P0: return None
    ka = int(np.argmax(dp > 0.05*dp.max()))
    k = int(np.argmax(dp)); ta = t[ka]
    w = (t >= ta) & (t <= ta + 1.5e-3)
    I15 = np.trapezoid(dp[w], t[w]) if t[-1] > ta + 1.5e-3 else None
    # humps: local maxima of the 0.05 ms smoothed trace 0.2-2.5 ms after arrival, rising > 3% of the peak
    m = (t > ta + 0.2e-3) & (t < ta + 2.5e-3)
    humps = 0
    if m.sum() > 20:
        tt, pp = t[m], dp[m]
        g = np.linspace(tt[0], tt[-1], 500); s = np.interp(g, tt, pp)
        s = np.convolve(s, np.ones(10)/10, "same")
        i = 1
        while i < len(s) - 1:
            if s[i] > s[i-1] and s[i] >= s[i+1]:
                lo = s[max(0, i-50):i].min()
                if s[i] - lo > 0.03*dp.max(): humps += 1
            i += 1
    return dict(ta=ta, pk=dp.max(), I15=I15, humps=humps, done=t[-1] > ta + 2.5e-3)

res = {}
for c in ("ctr", "top"):
    log = os.path.join(D, c, "log.blastEulerFoam")
    t = float(tail_val(log, "^Time = ").split()[-1]); dt = float(tail_val(log, "^deltaT = ").split()[-1])
    cells = tail_val(log, "(Refined|Unrefined) from").split()[-2] if tail_val(log, "(Refined|Unrefined) from") else "?"
    end = float(re.search(r"^endTime\s+([\d.eE+-]+);", open(os.path.join(D, c, "system/controlDict")).read(), re.M).group(1))
    err = subprocess.run(["grep", "-a", "-c", "-E", r"FATAL|Floating point|\bnan\b|NaN", log], capture_output=True, text=True).stdout.strip()
    done = subprocess.run(["grep", "-a", "-q", "^End", log]).returncode == 0
    alive = subprocess.run(["pgrep", "-f", f"blastEulerFoam.*"], capture_output=True, text=True).stdout.strip() != ""
    age = now - os.path.getmtime(log)
    prev = st.get(c)
    rate = (t - prev["t"])/((now - prev["now"])/3600) if prev and now > prev["now"] else None
    eta = (end - t)/rate if rate and rate > 0 else None
    st[c] = {"t": t, "now": now}
    f = []
    if err not in ("0", ""): f.append("ERROR")
    if age > 1800 and not done: f.append("LOG-STALE")
    if dt < 1e-9: f.append("DTLOW")
    if rate is not None and rate <= 0 and not done: f.append("STALL")
    if eta and eta > 72: f.append("SLOW(ETA>72h)")
    flags += [f"{c}:{x}" for x in f]
    out.append(f"{c}: t {t*1e3:.3f}/{end*1e3:.1f} ms, dt {dt:.2e}, cells {cells}, "
               + (f"rate {rate*1e3:.3f} ms/h, ETA {eta:.1f} h" if eta else "rate n/a")
               + (" DONE" if done else "") + (" " + " ".join(f) if f else ""))
    pos, d = probes(c)
    for j, (x, y) in enumerate(pos):
        r = round((x*x + y*y)**0.5, 3); ang = round(np.degrees(np.arctan2(y, x)), 1)
        fe = feats(d[:, 0], d[:, 1 + j])
        if fe: res[(c, ang, r)] = fe
# physics
angles = sorted({a for (_, a, _) in res}, reverse=True)
radii = sorted({r for (_, _, r) in res})
for r in radii:
    cs = [res[("ctr", a, r)] for a in angles if ("ctr", a, r) in res]
    if len(cs) == len(angles):
        tas = np.array([x["ta"] for x in cs]); pks = np.array([x["pk"] for x in cs])
        spread_ta = (tas.max() - tas.min())/tas.mean(); spread_pk = (pks.max() - pks.min())/pks.mean()
        dev1d = tas.mean()/TA1D.get(r, np.nan) - 1
        line = f"ctr r {r:g} m: ta {tas.mean()*1e3:.3f} ms (1D {dev1d:+.1%}), spread over angles ta {spread_ta:.1%} pk {spread_pk:.1%}, pk {pks.mean()/1e3:.1f} kPa"
        if spread_ta > 0.02 or spread_pk > 0.10: flags.append(f"ctr-ANISO r{r:g}")
        if abs(dev1d) > 0.05: flags.append(f"ctr-vs-1D r{r:g} {dev1d:+.1%}")
        out.append(line)
    for a in angles:
        if ("ctr", a, r) in res and ("top", a, r) in res:
            A, B = res[("ctr", a, r)], res[("top", a, r)]
            s = (f"  top-ctr ang {a:+6.1f} r {r:g}: dta {1e6*(B['ta']-A['ta']):+6.1f} us, pk {B['pk']/A['pk']-1:+.1%}"
                 + (f", I1.5ms {B['I15']/A['I15']-1:+.1%}" if A['I15'] and B['I15'] else "")
                 + (f", humps top {B['humps']} ctr {A['humps']}" if A['done'] and B['done'] else ""))
            out.append(s)
json.dump(st, open(STATE, "w"))
print(time.strftime("%a %H:%M"), "| flags:", ", ".join(flags) if flags else "none")
print("\n".join(out))
