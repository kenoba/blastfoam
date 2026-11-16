#!/usr/bin/env python3
"""Afterburning budget of an FOI 2013 case: gas species masses, chemical
energy released (formation enthalpies at 298.15 K from thermo.nasa9Gas) and
the fraction of the combustion potential (C -> CO2, H -> H2O) used.
usage: afterburn.py <case> <time> [<time> ...]  (fields reconstructed, Vc written)"""
import sys, os, re, math
import numpy as np
R = 8.314462618
case = sys.argv[1]; times = sys.argv[2:]

def field(t, name):
    s = open(os.path.join(case, t, name)).read()
    m = re.search(r"internalField\s+uniform\s+([^;]+);", s)
    if m: return float(m.group(1))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(", s)
    n = int(m.group(1)); body = s[m.end():]
    return np.array(body[:body.index(")")].split(), dtype=float)[:n]

th = open(os.path.join(case, "constant", "thermo.nasa9Gas")).read()
def hf(sp):
    m = re.search(r"\n\"?" + re.escape(sp) + r"\"?\s*\{.*?molWeight\s+([\d.eE+-]+);.*?coeffs\s*\(\s*\(([^)]*)\)", th, re.S)
    W = float(m.group(1)); a = [float(x) for x in m.group(2).split()]; T = 298.15
    h = (-a[0]/T**2 + a[1]*math.log(T)/T + a[2] + a[3]*T/2 + a[4]*T**2/3 + a[5]*T**3/4 + a[6]*T**4/5 + a[7]/T)*R*T
    return h/W*1e3, W   # J/kg, g/mol
sp = ["CO","C(gr)","CO2","H2","H2O","O2","N2","OH","H","O","NO","CH4","NH3","CH2O","HCO","HO2","H2O2","C2H6","CH3","CH3O","N"]
atoms = {"CO":(1,0),"C(gr)":(1,0),"CO2":(1,0),"H2":(0,2),"H2O":(0,2),"OH":(0,1),"H":(0,1),"CH4":(1,4),"NH3":(0,3),
         "CH2O":(1,2),"HCO":(1,1),"HO2":(0,1),"H2O2":(0,2),"C2H6":(2,6),"CH3":(1,3),"CH3O":(1,3)}
H = {s: hf(s) for s in sp}
hCO2 = H["CO2"][0]*44.0095e-3; hH2O = H["H2O"][0]*18.0153e-3   # J/mol
rows = []
for t in times:
    V = field(t, "Vc"); ar = field(t, "alphaRho.gas")
    m = {s: float(np.sum(ar*field(t, s+".gas")*V)) for s in sp}
    mc = float(np.sum(field(t, "alphaRho.compBAl")*V))
    Hchem = sum(m[s]*H[s][0] for s in sp)
    # enthalpy after complete combustion of C and H (O2 from the air), N as is
    nC = sum(m[s]/H[s][1]*1e3*atoms.get(s,(0,0))[0] for s in sp)
    nH = sum(m[s]/H[s][1]*1e3*atoms.get(s,(0,0))[1] for s in sp)
    Hfuel = sum(m[s]*H[s][0] for s in sp if s in atoms)
    pot = Hfuel - (nC*hCO2 + nH/2*hH2O)
    rows.append((t, m, mc, Hchem, pot, nC))
t0 = rows[0]
print(f"case {case}")
print(f"{'t [ms]':>7} {'compBAl':>9} {'CO':>9} {'C(gr)':>9} {'H2':>9} {'CO2':>9} {'H2O':>9} {'O2':>9} {'released':>10} {'potential':>10} {'used':>6}   (kg, J; per domain)")
for t, m, mc, Hc, pot, nC in rows:
    rel = t0[3] - Hc
    print(f"{float(t)*1e3:7.2f} {mc:9.3g} {m['CO']:9.3g} {m['C(gr)']:9.3g} {m['H2']:9.3g} {m['CO2']:9.3g} {m['H2O']:9.3g} {m['O2']:9.3g} {rel:10.4g} {pot:10.4g} {rel/t0[4]*100:5.1f}%")
print(f"carbon [mol]: {t0[5]:.5g} at {t0[0]}, {rows[-1][5]:.5g} at {rows[-1][0]}")
