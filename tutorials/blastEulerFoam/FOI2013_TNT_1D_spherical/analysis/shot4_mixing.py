#!/usr/bin/env python3
"""FOI shot 4 (TNT): effect of the subgrid mixing coefficient C.
(a) secondary shock after the main shock at 2.83 and 4.47 m (1D, 7 m mesh)
    against C, with the measured values at P1/P2;
(b) main shock at 4.47 m against C: arrival shift and positive impulse
    relative to C = 0.005;
(c, d) ground overpressure at P2/P3 (2D): measurement, laminar, C = 0.005
    and C = 0.0125.
    shot4_mixing.py <out.pdf>    (run directories from FOI2013_RUNS)
Secondary shock: the largest pressure rise in 30 us more than 0.8 ms after
arrival; main shock: first 2 kPa rise, positive impulse to the first zero
crossing."""
import os, re, sys
import numpy as np
sys.path.insert(0, os.environ.get("BLAST_DIR", "/opt/OpenFOAM/blastfoam")
                + "/tutorials/blastEulerFoam/FOI2013_TNT_1D_spherical")
import extract_shots as ex
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
trap = getattr(np, "trapezoid", None) or np.trapz
R = os.environ.get("FOI2013_RUNS", os.getcwd())
ONED = [(0.0, "secshock/s4_lam"), (0.005, "alburn/s4_base"), (0.0075, "secshock/s4_c0.0075"),
        (0.01, "secshock/s4_c0.01"), (0.015, "secshock/s4_c0.015"), (0.03, "secshock/s4_c0.03")]
TWOD = [("laminar", "peak/L3", "C2", "--"), ("C = 0.005", "al2d/s4_mix_L3", "C0", "-."),
        ("C = 0.0125", "c0125/s4_2D", "C3", "-")]
MEAS_SEC = {2.83: 2.66, 4.47: 3.01}   # measured at P1/P2 (2/4 m on the ground)


def probes(path):
    xs = [tuple(float(v) for v in re.search(r"\(([^)]*)\)", l).group(1).split())
          for l in open(path) if l.startswith("# Probe")]
    d = np.loadtxt(path, comments="#")
    return d[:, 0], d[:, 1:], xs


def col(xs, x, y=None):
    return min(range(len(xs)), key=lambda i: abs(xs[i][0] - x) + (abs(xs[i][1] - y) if y is not None else 0))


def metrics(t, p):
    y = p - p[0]
    a = np.argmax(y > 2000); ta = t[a]; k = a
    while k < len(y) and y[k] > 0:
        k += 1
    tg = np.arange(ta + 0.8e-3, min(t[-1], ta + 6e-3), 5e-6)
    yg = np.interp(tg, t, y); r = yg[6:] - yg[:-6]
    return ta, trap(y[a:k], t[a:k]), tg[np.argmax(r)] - ta


def main():
    out = sys.argv[1]
    fig, ax = plt.subplots(2, 2, figsize=(10, 7.5))
    res = {}
    for C, c in ONED:
        t, p, xs = probes(f"{R}/{c}/postProcessing/pressureProbes/0/p")
        res[C] = {r: metrics(t, p[:, col(xs, r)]) for r in (2.828, 4.472)}
    Cs = [C for C, _ in ONED]
    a = ax[0, 0]
    for r, m, lab in ((2.828, "o", "2.83 m"), (4.472, "s", "4.47 m")):
        a.plot(Cs, [res[C][r][2]*1e3 for C in Cs], m + "-", label=f"model, {lab}")
        a.axhline(MEAS_SEC[round(r, 2)], color=a.lines[-1].get_color(), ls=":", lw=1,
                  label=f"measured, P{1 if r < 3 else 2}")
    a.axvline(0.0125, color="0.6", lw=0.8)
    a.set_xlabel("C"); a.set_ylabel("secondary shock after main shock [ms]")
    a.set_title("(a) secondary shock (1D)", fontsize=10); a.legend(fontsize=8)
    a = ax[0, 1]
    ref = res[0.005][4.472]
    a.plot(Cs, [(res[C][4.472][0] - ref[0])*1e6 for C in Cs], "o-", label="arrival shift [$\\mu$s]")
    a2 = a.twinx()
    a2.plot(Cs, [100*(res[C][4.472][1]/ref[1] - 1) for C in Cs], "s-", color="C1", label="impulse [%]")
    a.set_xlabel("C"); a.set_ylabel("arrival shift at 4.47 m [$\\mu$s]")
    a2.set_ylabel("positive impulse at 4.47 m [%]", color="C1")
    a.set_title("(b) main shock relative to C = 0.005 (1D)", fontsize=10)
    a.axvline(0.0125, color="0.6", lw=0.8)
    tm, pm, _ = ex.load(4)
    for k, (r, win, ylim) in enumerate(((4.0, (7.0, 13.0), (-15, 110)), (6.0, (11.8, 17.5), (-10, 50)))):
        a = ax[1, k]
        a.plot(tm*1e3, np.convolve(pm[:, k + 1], np.ones(20)/20, "same")/1e3, color="k", lw=0.8,
               label="measured")
        for lab, c, colr, ls in TWOD:
            t, p, xs = probes(f"{R}/{c}/postProcessing/pressureProbes/0.0015/p")
            a.plot(t*1e3, (p[:, col(xs, r, 0.001)] - 101325)/1e3, color=colr, ls=ls, lw=1.1, label=lab)
        a.set_xlim(*win); a.set_ylim(*ylim); a.axhline(0, color="0.8", lw=0.5)
        a.set_xlabel("time after detonation [ms]"); a.set_ylabel("overpressure [kPa]")
        a.set_title(f"({'cd'[k]}) P{k + 2} ({r:.0f} m), 2D", fontsize=10)
    ax[1, 0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(out)
    for C in Cs:
        print(f"C {C:6.4f}: " + "  ".join(f"{r} m: arr {res[C][r][0]*1e3:.3f} ms I {res[C][r][1]:.1f} Pa s sec +{res[C][r][2]*1e3:.2f} ms" for r in (2.828, 4.472)))


if __name__ == "__main__":
    main()
