#!/usr/bin/env python3
"""FOI 2013 shot 4 1D case with subgrid turbulent mixing of products and air
(fvModel turbulentMixingDiffusion in the gas species and energy).

D(t) = C * w dw/dt with the mixing-layer width of Balakrishnan et al. (2010,
Shock Waves 20:147, 3D TNT): w/r0 = max(7.4 tau, 19.7 tau^1.56),
tau = t[ms]/W^(1/3), up to tau1 (end of their implosion phase, 0.9);
afterwards D decays as (tau1/tau)^late.

usage: setup_mix.py --out DIR --C 0.01 [--endTime 3e-3] [--R 4] [--tau1 0.9]
                    [--late 1] [--energy yes|no]
"""
import argparse
import math
import os
import re
import subprocess
import sys

BLAST = os.environ.get("BLAST_DIR", "/opt/OpenFOAM/blastfoam")
TUT = os.path.join(BLAST, "tutorials/blastEulerFoam/FOI2013_TNT_1D_spherical")

ap = argparse.ArgumentParser()
ap.add_argument("--out", required=True)
ap.add_argument("--C", type=float, required=True)
ap.add_argument("--endTime", type=float, default=3e-3)
ap.add_argument("--R", type=float, default=4.0)
ap.add_argument("--tau1", type=float, default=0.9)
ap.add_argument("--late", type=float, default=1.0)
ap.add_argument("--energy", default="yes")
ap.add_argument("--zMin", type=float, default=1e-3)
a = ap.parse_args()

subprocess.run([sys.executable, os.path.join(TUT, "setup_foi2013.py"), "--shot", "4",
                "--out", a.out, "--R", str(a.R), "--endTime", str(a.endTime)], check=True)

info = open(os.path.join(a.out, "charge.txt")).read().split()
W = float(info[info.index("m_total") + 1])
r0 = float(info[info.index("R") + 1])
W3 = W ** (1.0 / 3.0)


def width(tau):
    return r0 * max(7.4 * tau, 19.7 * tau ** 1.56)


def D(t):
    """C w dw/dt [m^2/s] at time t [s]"""
    if t <= 0:
        return 0.0
    tau = t * 1e3 / W3
    tc = min(tau, a.tau1)
    h = 1e-4 * tc
    dwdtau = (width(tc + h) - width(tc - h)) / (2 * h) if tc > h else 7.4 * r0
    dwdt = dwdtau * 1e3 / W3          # per second
    d = a.C * width(tc) * dwdt
    if tau > a.tau1:
        d *= (a.tau1 / tau) ** a.late
    return d


ts = sorted(set([0.0] + [a.endTime * (i / 400.0) ** 2 for i in range(1, 401)]))
table = "\n".join(f"            ({t:.6e} {D(t):.6e})" for t in ts)

# All species of the gas, the default specie (N2, the tracer) included: the
# gas receives the products by mass transfer, so its thermo solves every
# species and normalises; an undiffused species would be converted into the
# others (the model checks the list)
pp = open(os.path.join(a.out, "constant/physicalProperties.gas")).read()
spec = re.search(r"\nspecies\s*\((.*?)\);", pp, re.S).group(1).split()
spec = [s.strip('"') for s in spec]
YAir = 0.76
YProd = 0.18321973

fv = os.path.join(a.out, "constant/fvModels")
s = open(fv).read()
block = f"""
// Subgrid turbulent mixing of products and air in the gas phase
// (setup_mix.py --C {a.C:g} --tau1 {a.tau1:g} --late {a.late:g})
mixing
{{
    type            turbulentMixingDiffusion;
    phase           gas;
    species         ({' '.join('"' + x + '"' if '(' in x else x for x in spec)});
    energy          {a.energy};
    tracer          N2;
    YAir            {YAir};
    YProducts       {YProd};
    zMin            {a.zMin:g};
    D
    {{
        type            table;
        values
        (
{table}
        );
    }}
}}
"""
s = s.replace("// ************************************************************************* //",
              block + "\n// ************************************************************************* //", 1)
open(fv, "w").write(s)

# Tight species solves: the gas thermo normalises the mass fractions after the
# species are solved one by one, so the diffusion solver error would turn
# into element changes (turbulentMixingDiffusion warns above 1e-12)
fs = os.path.join(a.out, "system/fvSolution")
s, n = re.subn(r'("[^"\n]*\\.gas[^"\n]*"\s*\{[^}]*?tolerance\s+)[^;]+;', r"\g<1>1e-15;",
               open(fs).read())
if n != 1:
    sys.exit(f"{fs}: no single species solver entry to tighten")
open(fs, "w").write(s)

# Probes at Balakrishnan's Z = 0.80 and 1.38 m/kg^1/3 for this charge
cd = os.path.join(a.out, "system/controlDict")
s = open(cd).read()
extra = f"            ({0.80 * W3:.3f} 0 0)\n            ({1.38 * W3:.3f} 0 0)\n"
s = re.sub(r"(probeLocations\s*\(\s*\n)", r"\g<1>" + extra, s, count=1)
open(cd, "w").write(s)

print(f"W {W:.4f} kg, r0 {r0:.4f} m; D at 0.1/0.3/0.6/0.9/1.5/3 ms:",
      " ".join(f"{D(t * 1e-3):.3g}" for t in (0.1, 0.3, 0.6, 0.9, 1.5, 3.0)), "m2/s")
