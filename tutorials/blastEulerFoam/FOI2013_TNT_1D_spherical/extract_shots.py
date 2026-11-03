#!/usr/bin/env python3
"""Arrival, peak and positive impulse at P1-P3 (ground-flush, 2/4/6 m) and the
pyrometer maximum of the FOI June 2013 TNT and TNT/Al shots (Fedina,
FOI-R--3913--SE, Table 4.1; raw data (not in the repository) in
report/external_sources/data_sources/Katja_data/Measurements/TXT or $FOI2013_DATA, 1 MHz,
overpressure in Pa, pyrometer in deg C with a 650 C floor). Time zero is the
trigger-wire short. Shot N = 'Katja SkottN' (checked: 1.SET labels shot 1
'1kg TNT 0.15 m'; the TNT peaks match Table 4.2)."""
import os, sys
import numpy as np
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    os.environ.get("FOI2013_DATA", "../../../report/external_sources/data_sources/Katja_data/Measurements/TXT"))
SHOTS = {1: ('TNT', 0.9103, 0.15), 2: ('TNT', 0.9425, 0.5), 3: ('TNT', 0.9093, 1), 4: ('TNT', 0.9002, 2),
         5: ('TNT/Al', 1.0514, 2), 6: ('TNT/Al', 1.0146, 1), 7: ('TNT/Al', 1.0666, 0.5),
         8: ('TNT/Al', 1.0009, 0.15), 9: ('TNT', 0.913, 0.5), 10: ('TNT', 1.845, 1), 11: ('TNT', 1.87, 0.5)}


def load(n):
    lines = open(os.path.join(DATA, f"Katja Skott{n:03d}.txt"), encoding="latin1").read().split("\n")
    ih = [i for i, l in enumerate(lines[:20]) if l.startswith("Time,")][0]
    d = np.genfromtxt(lines[ih + 1:], delimiter=",", invalid_raise=False)
    tr = d[:, 5]
    k0 = int(np.argmax(np.abs(tr - np.median(tr[:1000])) > 1.0))
    t = d[:, 0] - d[k0, 0]
    p = d[:, 1:4] - np.median(d[:k0, 1:4], axis=0)
    return t, p, d[:, 4]


def metrics(t, p, thr):
    idx = np.where((t > 0.2e-3) & (p > thr))[0]
    ka = idx[0]
    while ka > 0 and p[ka - 1] > 0.25*thr:
        ka -= 1
    pk = p[(t >= t[ka]) & (t < t[ka] + 0.3e-3)].max()
    pm = np.convolve(p, np.ones(200)/200, "same")      # 0.2 ms running mean
    j = ka + 50 + int(np.argmax(pm[ka + 50:] < 0))
    return t[ka], pk, np.trapezoid(p[ka:j], t[ka:j]), t[j] - t[ka]


if __name__ == "__main__":
    print("shot charge  m_kg   HOB |  P1 ta[ms] pk[kPa] I+[Pa s] td[ms] |  P2 ... |  P3 ... | pyro max [C] t[ms]")
    for n, (ch, m, h) in SHOTS.items():
        t, p, py = load(n)
        s = f"{n:4d} {ch:<7}{m:6.3f} {h:5.2f} |"
        for c, thr in enumerate((20e3, 8e3, 4e3)):
            ta, pk, I, td = metrics(t, p[:, c], thr)
            s += f" {ta*1e3:7.3f} {pk/1e3:6.1f} {I:6.1f} {td*1e3:5.2f} |"
        k = int(np.argmax(py))
        print(s + f" {py[k]:6.0f} {t[k]*1e3:6.1f}")
