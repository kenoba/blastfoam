#!/usr/bin/env python3
"""Products isentrope of the FOI shot 4 model (eJWLC to the transfer pressure,
then the frozen Abel-Nobel/nasa9 gas) against the EuDetCode BKW isentrope
(extended below its last point with the same frozen gas).
    isentrope.py <case> <eudet dir>"""
import math, re, sys, os
import numpy as np
if not hasattr(np, "trapezoid"):  # numpy < 2.0
    np.trapezoid = np.trapz
case, eud = sys.argv[1:3]
R = 8.314462618
con = os.path.join(case, "constant")
txt = open(os.path.join(con, "physicalProperties.compBAl")).read()
prod = txt[txt.index("products"):]
g = lambda k: float(re.search(rf"\b{k}\s+([-\d.eE+]+);", prod).group(1))
A, B, C, R1, R2, om, rho0, E0, Tcj, rhocj = (g(k) for k in
    ("A", "B", "C", "R1", "R2", "omega", "rho0", "E0", "T_CJ", "rho_CJ"))
v0 = 1/rho0
ps = lambda V: A*np.exp(-R1*V) + B*np.exp(-R2*V) + C*V**(-1-om)
es = lambda V: v0*(A/R1*np.exp(-R1*V) + B/R2*np.exp(-R2*V) + C/om*V**(-om)) + E0
fv = open(os.path.join(con, "fvModels")).read()
pt = float(re.search(r"pThreshold\s+([\d.eE+-]+);", fv).group(1))
Y = {k.strip('"'): float(v) for k, v in re.findall(r"^\s+(\"?[\w()]+\"?)\s+([\d.eE+-]+);",
     fv[fv.index("targetFractions"):fv.index("isMoleFraction")], re.M)}
th = open(os.path.join(con, "thermo.nasa9Gas")).read()
coef = {}
for sp in Y:
    blk = th[th.index("\n" + ('"%s"' % sp if "(" in sp else sp) + "\n"):]
    blk = blk[:blk.index("\n}\n")]
    W = float(re.search(r"molWeight\s+([\d.eE+-]+);", blk).group(1))
    b = float(re.search(r"\bb\s+([\d.eE+-]+);", blk).group(1))
    Ts = [float(x) for x in re.search(r"Tswitch\s*\(([^)]*)\)", blk).group(1).split()]
    cs = [[float(x) for x in c.split()] for c in re.findall(r"\(([-\d.eE+ ]+)\)", blk[blk.index("coeffs"):])]
    coef[sp] = (W, b, Ts, cs, re.search(r"condensed\s+(yes|true|on)", blk) is not None)
def h(sp, T):
    W, b, Ts, cs, cond = coef[sp]
    i = min(max(np.searchsorted(Ts, T) - 1, 0), len(cs) - 1); a = cs[i]
    return (-a[0]/T**2 + a[1]*math.log(T)/T + a[2] + a[3]*T/2 + a[4]*T**2/3
            + a[5]*T**3/4 + a[6]*T**4/5 + a[7]/T)*R*T/W*1e3
Rg = sum(Y[s]*R/coef[s][0]*1e3 for s in Y if not coef[s][4])
bm = sum(Y[s]*coef[s][1] for s in Y)
eg = lambda T: sum(Y[s]*(h(s, T) - h(s, 298.15)) for s in Y) - Rg*T
pg = lambda rho, T: rho*Rg*T/(1 - bm*rho)
def Tof(e, T=2000.):
    for _ in range(60):
        d = (eg(T) - e)/((eg(T + .5) - eg(T - .5))); T -= d
        if abs(d) < 1e-8: break
    return T
def gas_isentrope(v, T, pend, n=4000):
    """frozen gas: de = -p dv, from (v, T) to pend; returns arrays v, p, e, T"""
    out = [(v, pg(1/v, T), eg(T), T)]
    f = math.exp(math.log(10)/200)
    while out[-1][1] > pend:
        v1 = v*f; T1 = T
        for _ in range(30):  # trapezoid in p
            T1n = Tof(eg(T) - 0.5*(pg(1/v, T) + pg(1/v1, T1))*(v1 - v), T1)
            if abs(T1n - T1) < 1e-9: break
            T1 = T1n
        v, T = v1, T1n; out.append((v, pg(1/v, T), eg(T), T))
    return np.array(out)
# model: eJWLC from V_CJ to transfer
Vcj = rho0/rhocj
Vt = Vcj
while ps(Vt) > pt: Vt *= 1.0005
Vs = np.geomspace(Vcj, Vt, 3000)
jw = np.column_stack([Vs*v0, ps(Vs), es(Vs)])
et = es(Vt); Tt = Tof(et)
print(f"eJWLC CJ: p {ps(Vcj)/1e9:.2f} GPa e_s {es(Vcj)/1e6:.4f} MJ/kg; transfer at p {pt/1e5:.0f} bar: "
      f"rho {rho0/Vt:.2f} kg/m3, e_s {et/1e6:.4f} MJ/kg")
print(f"frozen gas at the same (rho, e): T {Tt:.0f} K, p {pg(rho0/Vt, Tt)/1e5:.1f} bar (R {Rg:.1f}, b {bm:.3e})")
gm = gas_isentrope(Vt*v0, Tt, 101325)
model = np.vstack([jw, gm[1:, :3]])
# BKW
bkw = np.loadtxt(os.path.join(eud, "isentrope_output.dat"), comments="#")[:, [0, 1, 4, 2]]
vb, pb, eb, Tb = bkw[-1]
# extend with the frozen gas from the last BKW (v, p): T from the gas EOS at that (v, p)
Tx = pb*(1 - bm/vb)*vb/Rg
gx = gas_isentrope(vb, Tx, 101325)
print(f"BKW last point: p {pb/1e5:.1f} bar, T {Tb:.0f} K, rho {1/vb:.2f}; gas EOS T at that (rho,p) {Tx:.0f} K")
gx[:, 2] += eb - gx[0, 2]
bk = np.vstack([bkw[:, :3], gx[1:, :3]])
def at(arr, p):  # v and work from CJ at pressure p (log-interp)
    lp = np.log(arr[:, 1][::-1]); 
    v = np.exp(np.interp(math.log(p), lp, np.log(arr[:, 0][::-1])))
    e = np.interp(math.log(p), lp, arr[:, 2][::-1])
    return v, arr[0, 2] - e
print(f"{'p [bar]':>9} {'v BKW':>9} {'v model':>9} {'W BKW':>8} {'W model':>8}  [MJ/kg from CJ]")
for p in (1e5*x for x in (10000, 1000, 300, 100, 30, 10, 3, 1.01325)):
    a, b = at(bk, p), at(model, p)
    print(f"{p/1e5:9.2f} {a[0]:9.5f} {b[0]:9.5f} {a[1]/1e6:8.3f} {b[1]/1e6:8.3f}  dW {100*(b[1]/a[1]-1):+5.1f}%")
# work to the ambient state including -p0 dv (useful blast work ~ int (p-p0) dv)
for name, arr in (("BKW", bk), ("model", model)):
    v, p = arr[:, 0], arr[:, 1]
    print(f"{name}: int (p - p0) dv from CJ to p0 = {np.trapezoid(p - 101325, v)/1e6:.3f} MJ/kg, "
          f"v(p0) {v[-1]:.4f} m3/kg")
print(f"model gas T at 1 atm {gm[-1,3]:.0f} K, BKW-extension T at 1 atm {gx[-1,3]:.0f} K")
