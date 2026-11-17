import os, re, sys, numpy as np
sys.path.insert(0, os.environ.get("BLAST_DIR", "/opt/OpenFOAM/blastfoam") + "/tests/1D_spherical")
import kingery_bulmash as kb
def rd(f):
    s = open(f).read(); m = re.search(r"nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    n = int(m.group(1)); return np.array(s[m.end():s.index(")", m.end())].split(), float)[:n]
x, y, p = rd("0.018/Ccx"), rd("0.018/Ccy"), rd("0.018/pTemporalMax") - 101325
hob = 2.0
def near(px, py):
    i = np.argmin((x-px)**2 + (y-py)**2); return p[i], x[i], y[i]
W = 0.94
print("slant  free(2D,y=2)  KB(0.94kg)  1D-2mm | ground x  2D-ground  meas  2D/free  meas/free(1D)")
oned = {2.83: 80.3, 4.47: 35.2, 6.33: 20.3}
meas = {2: 440.9, 4: 99.9, 6: 45.1}
for g, R in ((2, 2.83), (4, 4.47), (6, 6.33)):
    pf = near(R, hob)[0]/1e3
    pk = kb.free_air_overpressure(R/W**(1/3))/1e3
    pg = near(g, 0.001)[0]/1e3
    print(f"{R:5.2f} {pf:10.1f} {pk:10.1f} {oned[R]:9.1f} | {g:6.1f} {pg:9.1f} {meas[g]:7.1f} {pg/pf:7.2f} {meas[g]/oned[R]:7.2f}")
# ground profile and the free-field along y = hob
print("ground profile x: p [kPa]"); print(" ".join(f"{gx}:{near(gx,0.001)[0]/1e3:.0f}" for gx in np.arange(0.5, 7.01, 0.5)))
print("vertical above charge (r up to 3 m):", " ".join(f"{r}:{near(0.004, hob+r)[0]/1e3:.0f}" for r in (1, 2, 2.83, 3)))
print("horizontal y=hob:", " ".join(f"{r}:{near(r, hob)[0]/1e3:.0f}" for r in (1, 2, 2.83, 3, 4, 4.47, 5, 6, 6.33, 7)))
print("KB 0.94 kg:", " ".join(f"{r}:{kb.free_air_overpressure(r/W**(1/3))/1e3:.0f}" for r in (1, 2, 2.83, 3, 4, 4.47, 5, 6, 6.33, 7)))
