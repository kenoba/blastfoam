#!/usr/bin/env python3
"""
verify_condensed_species.py: condensed species in the Abel-Nobel gas
('condensed yes' in the equationOfState, e.g. C(gr) in detonation products).

Generates the nasa9 species file of the TNT products (FOI 2013 transfer
composition, 17.5% C(gr) by mass) with utilities/thermo_fitter/
nasa9_species.py, once with C(gr) condensed (default) and once as a gas
(--condensed-as-gas, the files before 2026-11-08), and runs the 0D
closed-volume isentropic expansion of applications/test/
Test-isentropicExpansion (built if missing): the compiled thermo classes
integrate de = -p dv from 100 kg/m^3 and 3000 K to 0.1 kg/m^3, the
temperature following from the energy with TRhoE as in the flow solver.

Independent reference (Python, from the NASA-9 coefficients of the file):
R_g = sum Y_k R_k over the gaseous species, b = sum Y_k b_k,
e(T) = sum Y_k (h_k(T) - h_k(298.15)) - R_g T, and the isentrope
cv(T) dT = -R_g T dv/(v - b), cv = de/dT (RK4 in ln v).

Pass criteria:
- species: condensed C(gr) has R_g = 0, Cv = Cp, p = 0 and Es = Hs; gas
  species R_g = R_u/W and Cp - Cv = R_g (1e-10);
- mixture: R_g, b, p(rho, T) and Es(T) as the reference (1e-10);
- expansion, both treatments: energy consistent with the temperature
  (|e - e_ref(T)| < 1e-8 of the energy change) and T on the reference
  isentrope within 1e-5;
- the condensed treatment gives a hotter isentrope than the gas treatment.
"""

import math
import os
import subprocess
import sys

ROOT = "/opt/OpenFOAM/blastfoam"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "condensed_species")
BASE = os.path.join(HERE, "../1D_spherical/ChargeCompBAl_1D_spherical_reacting/constant/thermo.compressibleGas")

# TNT products at transfer (FOI 2013 cases, fvModels targetFractions)
Y = {"H2": 0.00299992, "N2": 0.18321973, "CO": 0.27437665, "H2O": 0.12795739,
     "CO2": 0.21483229, "CH4": 0.02126888, "C(gr)": 0.17534515}
# the compiled mixture is normalised by the total mass fraction
Y = {k: v/sum(Y.values()) for k, v in Y.items()}
RHO0, T0, RHO1 = 100.0, 3000.0, 0.1
RHO_TEST, T_TEST = 10.0, 1500.0
# Universal gas constant of the compiled specie class (OpenFOAM's RR differs
# from CODATA 2018 in the 7th digit); set from the N2 output in check()
RU = None
TSTD = 298.15


def env(cmd):
    return ("source /opt/OpenFOAM/OpenFOAM-12/etc/bashrc >/dev/null 2>&1 && "
            "source /opt/OpenFOAM/blastfoam/etc/bashrc >/dev/null 2>&1 && " + cmd)


def parse_species(path):
    """name -> (W, b, condensed, Tswitch, coeffs) from a nasa9 species file."""
    import re
    txt = open(path).read()
    out = {}
    for name in Y:
        q = re.escape(f'"{name}"' if "(" in name else name)
        m = re.search(r"^" + q + r"\s*\{(.*?)^\}", txt, re.S | re.M)
        blk = m.group(1)
        W = float(re.search(r"molWeight\s+([^;]+);", blk).group(1))
        eos = re.search(r"equationOfState\s*\{([^}]*)\}", blk).group(1)
        b = float(re.search(r"\bb\s+([^;]+);", eos).group(1))
        cond = re.search(r"condensed\s+(yes|true|on)\s*;", eos) is not None
        Ts = [float(x) for x in re.search(r"Tswitch\s*\(([^)]*)\)", blk).group(1).split()]
        cblk = re.search(r"coeffs\s*\((.*?)\);", blk, re.S).group(1)
        co = [[float(x) for x in r.split()] for r in re.findall(r"\(([^()]*)\)", cblk)]
        out[name] = (W, b, cond, Ts, co)
    return out


class Mixture:
    def __init__(self, sp):
        self.sp = sp
        self.R = {k: (0.0 if sp[k][2] else RU/sp[k][0]) for k in Y}
        self.Rg = sum(Y[k]*self.R[k] for k in Y)
        self.b = sum(Y[k]*sp[k][1] for k in Y)

    def _a(self, k, T):
        Ts, co = self.sp[k][3], self.sp[k][4]
        i = 0
        while i < len(co) - 1 and T > Ts[i + 1]:
            i += 1
        return co[i]

    def cp_k(self, k, T):
        a = self._a(k, T)
        return RU/self.sp[k][0]*(a[0]/T**2 + a[1]/T + a[2] + T*(a[3] + T*(a[4] + T*(a[5] + T*a[6]))))

    def ha_k(self, k, T):
        a = self._a(k, T)
        return RU/self.sp[k][0]*(-a[0]/T + a[1]*math.log(T) + a[7]
                                 + T*(a[2] + T*(a[3]/2 + T*(a[4]/3 + T*(a[5]/4 + T*a[6]/5)))))

    def hs_k(self, k, T):
        return self.ha_k(k, T) - self.ha_k(k, TSTD)

    def es_k(self, k, T):
        return self.hs_k(k, T) - self.R[k]*T

    def e(self, T):
        return sum(Y[k]*self.es_k(k, T) for k in Y)

    def cv(self, T):
        return sum(Y[k]*(self.cp_k(k, T) - self.R[k]) for k in Y)

    def p(self, rho, T):
        return rho*self.Rg*T/(1 - self.b*rho)

    def isentrope(self, rhos, n=200000):
        """T at the densities rhos (decreasing from RHO0), RK4 in x = ln v."""
        f = lambda x, T: -self.Rg*T*math.exp(x)/((math.exp(x) - self.b)*self.cv(T))
        x, T = math.log(1/RHO0), T0
        out = []
        xs = [math.log(1/r) for r in rhos]
        h = (math.log(1/RHO1) - x)/n
        j = 0
        while j < len(xs) and xs[j] <= x + 1e-15:
            out.append(T)
            j += 1
        for _ in range(n):
            k1 = f(x, T)
            k2 = f(x + h/2, T + h/2*k1)
            k3 = f(x + h/2, T + h/2*k2)
            k4 = f(x + h, T + h*k3)
            Tn = T + h/6*(k1 + 2*k2 + 2*k3 + k4)
            while j < len(xs) and xs[j] <= x + h + 1e-12:
                # linear in x within the step
                w = (xs[j] - x)/h
                out.append(T + w*(Tn - T))
                j += 1
            x, T = x + h, Tn
        while j < len(xs):
            out.append(T)
            j += 1
        return out


def run(species_file):
    mix = "(" + " ".join(f'({k if "(" not in k else chr(34) + k + chr(34)} {v})' for k, v in Y.items()) + ")"
    cmd = (f"Test-isentropicExpansion -file {species_file} -mixture '{mix}' "
           f"-rho0 {RHO0} -T0 {T0} -rho1 {RHO1} -nSteps 20000")
    out = subprocess.run(["bash", "-c", env(cmd)], capture_output=True, text=True, cwd=WORK)
    if out.returncode != 0:
        print(out.stdout[-2000:], out.stderr[-2000:])
        sys.exit("FAIL: Test-isentropicExpansion failed")
    specie, mixture, states = {}, None, []
    for l in out.stdout.splitlines():
        r = l.split()
        if not r:
            continue
        if r[0] == "specie":
            specie[r[1].strip('"')] = [float(x) for x in r[2:]]
        elif r[0] == "mixture":
            mixture = [float(x) for x in r[1:]]
        elif r[0] == "state":
            states.append([float(x) for x in r[1:]])
    return specie, mixture, states


def check(label, species_file, condensed):
    global RU
    sp = parse_species(species_file)
    specie, mixture, states = run(species_file)
    RU = specie["N2"][0]*sp["N2"][0]
    m = Mixture(sp)
    ok = True
    msgs = []

    def rel(a, b, scale=None):
        return abs(a - b)/(scale if scale else max(abs(b), 1e-300))

    # species at the test state: Rg b Cp Cv Es p
    worst = 0.0
    for k in Y:
        Rg, b, Cp, Cv, Es, p = specie[k]
        Rk = m.R[k]
        errs = [abs(Rg - Rk)/(RU/sp[k][0]), rel(b, sp[k][1]), rel(Cp, m.cp_k(k, T_TEST)),
                abs((Cp - Cv) - Rk)/Cp, abs(Es - m.es_k(k, T_TEST))/abs(m.hs_k(k, T_TEST)),
                abs(p - RHO_TEST*Rk*T_TEST/(1 - sp[k][1]*RHO_TEST))/(RHO_TEST*RU/sp[k][0]*T_TEST)]
        worst = max(worst, max(errs))
        if k == "C(gr)":
            msgs.append(f"  C(gr): Rg {Rg:.6g}, Cp - Cv {Cp - Cv:.3g}, p {p:.6g} Pa at "
                        f"{RHO_TEST:g} kg/m3 and {T_TEST:g} K")
            if condensed and not (Rg == 0 and Cp == Cv and p == 0):
                ok = False
                msgs.append("  FAIL: C(gr) not condensed")
    msgs.append(f"  species Rg, b, Cp, Cp - Cv, Es, p vs reference: max rel. diff {worst:.1e}")
    ok = ok and worst < 1e-10

    Rg, b, Cp, Cv, Es, p = mixture
    errs = [rel(Rg, m.Rg), rel(b, m.b), rel(Cv, m.cv(T_TEST)),
            abs(Es - m.e(T_TEST))/abs(m.e(T_TEST) - m.e(TSTD)), rel(p, m.p(RHO_TEST, T_TEST))]
    msgs.append(f"  mixture Rg {Rg:.6g} J/kg/K, b {b:.6g} m3/kg; Rg, b, Cv, Es, p vs reference: "
                f"max rel. diff {max(errs):.1e}")
    ok = ok and max(errs) < 1e-10

    rhos = [s[0] for s in states]
    Tref = m.isentrope(rhos)
    de = abs(states[-1][3] - states[0][3])
    eErr = max(abs(s[3] - m.e(s[1]))/de for s in states)
    TErr = max(abs(s[1] - t)/t for s, t in zip(states, Tref))
    pErr = max(abs(s[2] - m.p(s[0], s[1]))/s[2] for s in states)
    msgs.append(f"  expansion {RHO0:g} -> {RHO1:g} kg/m3 from {T0:g} K: |e - e(T)|/de {eErr:.1e}, "
                f"T vs isentrope {TErr:.1e}, p vs EOS {pErr:.1e}")
    ok = ok and eErr < 1e-8 and TErr < 1e-5 and pErr < 1e-10

    def T_at(r):
        for s in states:
            if s[0] <= r*(1 + 1e-9):
                return s[1]
        return states[-1][1]
    Ts = {r: T_at(r) for r in (10, 1, RHO1)}
    msgs.append("  T at " + ", ".join(f"{r:g} kg/m3: {t:.0f} K" for r, t in Ts.items()))
    print(f"{label}:")
    print("\n".join(msgs))
    return ok, Ts


def main():
    os.makedirs(WORK, exist_ok=True)
    r = subprocess.run(["bash", "-c", env("which Test-isentropicExpansion")], capture_output=True)
    if r.returncode != 0:
        subprocess.run(["bash", "-c", env(f"wmake {ROOT}/applications/test/Test-isentropicExpansion")],
                       check=True, capture_output=True)
    gen = os.path.join(ROOT, "utilities/thermo_fitter/nasa9_species.py")
    files = {}
    for label, extra in (("condensed", []), ("gas", ["--condensed-as-gas"])):
        files[label] = os.path.join(WORK, f"thermo.nasa9Gas.{label}")
        subprocess.run([sys.executable, gen, "--base", BASE, "--species", *Y, "-o", files[label], *extra],
                       check=True, capture_output=True)
    ok_c, Tc = check("C(gr) condensed", files["condensed"], True)
    ok_g, Tg = check("C(gr) as a gas", files["gas"], False)
    hotter = all(Tc[r] > Tg[r] for r in Tc)
    print(f"condensed isentrope hotter than the gas one at all densities: {hotter}")
    ok = ok_c and ok_g and hotter
    print("SUCCESS: condensed species in the Abel-Nobel gas" if ok
          else "FAILURE: condensed species in the Abel-Nobel gas")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
