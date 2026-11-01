#!/usr/bin/env python3
"""Peak overpressure shape factor of centre-initiated cylindrical charges in
free air, Sp = p_cylinder/p_sphere of equal mass, from the AUTODYN fit of
Gao, Kong, Fang, Hong & Wang, Defence Technology 18 (2022) 1662-1678,
Eq. (9) with Table 2 (TNT, L/D >= 2, 0.3 m/kg^1/3 <= Z <= Zp, Eq. (7)).
The sphere is the Henrych (1979) TNT fit, which Gao et al. found to match
their spherical runs (their Fig. 6).

Angles here follow the G3a probes: 0 deg side-on (mid-plane), 90 deg
end-on (axis); Gao et al. use the opposite (0 axial, 90 radial).

Usage: cylinder_blast_reference.py <g3_early_table.txt> [TNT equivalence]
compares the L/D >= 2 rows of the table written by g3_early.py
(results/G3a_cylinders_early/) with the fit. Z uses the charge mass times
the TNT equivalence (default 1; 0.953 kg Comp-B/Al).
"""
import sys

import numpy as np

T = {  # range -> i -> [a_i1, a_i2, (a_ij1, a_ij2) j=1..4]
 0: {1: [10150, -17744, 1804, 3174, -18972, 24302, 14156, -16584, -10943, 14167],
     2: [-4209, 12410, -516, -4237, 9042, -14271, -8585, 10029, 7534, -10020],
     3: [781, -2006, -1424, 2004, 463, 525, 88, -360, -869, 1233]},
 1: {1: [-448, 396, -1105, 776, 2343, -1783, -1622, 1317, 828, -631],
     2: [1247, -1198, 6717, -3749, -10495, 7386, 7820, -5600, -3600, 2552],
     3: [1143, 93, -5571, 3299, 6144, -4862, -5569, 3920, 2295, -1675]},
 2: {1: [-59, 169, -555, 449, 945, -886, -626, 561, 272, -252],
     2: [1897, -1431, 6481, -5095, -11258, 9292, 7652, -6036, -3419, 2704],
     3: [-4387, 3190, -8614, 10132, 19388, -18959, -12744, 12005, 5995, -5410]}}
def Z1(ld): return 0.613 + 5.411*np.exp(-0.962*ld)
def Z2(ld): return 10.116 + 49.5*np.exp(-0.341*ld) - 51.942*np.exp(-0.253*ld)
def Zp(ld): return 4.749 + 22.033*np.exp(-0.758*ld)
def p(Z, theta, ld):
    r = 0 if Z <= Z1(ld) else (1 if Z <= min(Z2(ld), Zp(ld)) else 2)
    L, th, s = np.log(ld), np.radians(theta), 0.0
    for i, a in T[r].items():
        c = a[0] + a[1]*L + sum((a[2*j] + a[2*j+1]*L)*np.sin(j*th) for j in range(1, 5))
        s += c/Z**i
    return s
def henrych(Z):  # kPa, sphere TNT
    if Z < 0.3: b = 14.072/Z + 5.540/Z**2 - 0.357/Z**3 + 0.00625/Z**4
    elif Z <= 1: b = 6.194/Z - 0.326/Z**2 + 2.132/Z**3
    else: b = 0.662/Z + 4.05/Z**2 + 3.288/Z**3
    return 100*b


def shape_factor(Z, angle, ld):
    """Sp at scaled distance Z [m/kg^1/3] and probe angle [deg, 0 side-on]."""
    return p(Z, 90.0 - angle, ld)/henrych(Z)


if __name__ == "__main__":
    table = sys.argv[1]
    W = 0.953*(float(sys.argv[2]) if len(sys.argv) > 2 else 1.0)
    print("L/D  r[m]  angle    Z    Sp(run)  Sp(Gao)")
    for line in open(table):
        f = line.split()
        if len(f) != 7 or not f[0].isdigit():
            continue
        ld, r, a, sp = int(f[0]), float(f[1]), int(f[2]), float(f[5])
        Z = r/W**(1/3)
        if ld < 2 or Z < 0.3 or Z > Zp(ld):
            continue
        print(f"{ld:3d} {r:5.2f} {a:5d} {Z:6.3f} {sp:8.2f} {shape_factor(Z, a, ld):8.2f}")
