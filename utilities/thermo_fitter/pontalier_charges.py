#!/usr/bin/env python3
"""
pontalier_charges.py: charges of the Pontalier et al. (2021) free-field
trials (Shock Waves 31:203-230, doi:10.1007/s00193-021-00994-z) and their
EuDetCode products-only fits.

Charge (Sect. 2, Table 1): gelled nitromethane (95.6 NM / 4 PMMA / 0.4 glass
microballoons by mass) with embedded particles (Al H-10/H-50/H-95, Mg, glass
beads) in a 12.7 cm glass bulb, and a central C-4 booster of 86-93 g.

Model conventions (the same as the Comp-B/Al cases):
- the booster is lumped into the charge composition (C-4 as in the EuDetCode
  formulation: RDX/DOA/PIB/mineral oil 91/5.3/2.1/1.6);
- the microballoons are voids: their volume counts in the density, their
  mass (0.4% of the gel) is left out;
- Al: 90% of the Al is inert in the fit (--inert "Al(s):0.9") and forms the
  separate Al phase; the rest is in the products;
- glass beads: fitted as fully inert Al (--inert "Al(s):1.0"), so the fit
  has the products of the matrix alone with the beads in thermal
  equilibrium at CJ; in the solver the particle phase gets the glass density
  and heat capacity and does not ignite.
PMMA is not in the EuDetCode database: it is added to a copy of Data/ in a
temporary run directory (C5H8O2, h0 = -441 kJ/mol: liquid MMA -386 plus
polymerisation -55).

    python3 pontalier_charges.py [--run TRIAL ...] [--eos BKW]

Without --run it prints the charge table. With --run it runs CJDetonation
for the trials into eudet_pontalier/<trial>/ next to this script.
"""

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
BLAST = os.environ.get("BLAST_DIR", os.path.normpath(os.path.join(HERE, "../..")))
EUDET = os.path.join(BLAST, "external", "eudetcode")

# Densities [kg/m3]
RHO = {"NM": 1137.0, "PMMA": 1180.0, "GMB": 125.0, "C4": 1590.0,
       "Al": 2700.0, "glass": 2500.0, "Mg": 1740.0}
GEL = {"NM": 0.956, "PMMA": 0.040, "GMB": 0.004}
C4 = {"RDX": 0.91, "DOA": 0.053, "PIB": 0.021, "Mineral-Oil": 0.016}

# Table 1: trial -> (material, wt%, d [um], booster [g], mixture [g], p0 [kPa], T0 [C])
TRIALS = {
    "0": ("baseline", 0, 0, 87.2, 979.2, 92.79, 24),
    "1": ("Al", 10, 13, 88.9, 1039.8, 93.12, 26),
    "2": ("Al", 30, 13, 92.2, 1183.5, 92.71, 26),
    "3": ("Al", 40, 13, 91.3, 1257.0, 93.59, 26),
    "4": ("Al", 10, 54, 87.6, 1035.3, 92.96, 25),
    "5": ("Al", 20, 54, 90.8, 1106.4, 92.79, 24),
    "6": ("Al", 30, 54, 88.89, 1179.5, 92.89, 25),
    "7": ("Al", 40, 54, 87.7, 1281.3, 93.67, 22),
    "8": ("Al", 40, 54, 86.0, 1303.4, 92.96, 25),
    "9": ("Al", 10, 114, 90.0, 1044.6, 93.19, 25),
    "10": ("Al", 30, 114, 85.7, 1175.4, 92.76, 26),
    "11": ("Al", 40, 114, 92.2, 1298.5, 93.70, 25),
    "12": ("Al", 50, 114, 87.7, 1407.0, 93.20, 26),
    "17": ("glass", 10, 66, 92.7, 1038.5, 93.21, 23),
    "18": ("glass", 30, 66, 87.4, 1169.9, 92.75, 22),
    "19": ("glass", 50, 66, 87.3, 1378.5, 93.19, 22),
}

PMMA_ENTRY = {
    "name": "PMMA",
    "description": "Poly(methyl methacrylate), monomer unit C5H8O2 "
                   "(pontalier_charges.py): h0 = Hf(MMA, l) -386 kJ/mol "
                   "+ polymerisation -55 kJ/mol.",
    "aliases": ["Polymethylmethacrylate"],
    "elements": {"C": 5, "H": 8, "O": 2},
    "h0": -441000.0,
    "constCpParams": {"Cp": 142.0},
    "v0": 100.12e-3/1180.0,
    "murnaghanParams": {"B0": 5.0e9, "n": 6.0},
}


def charge(trial):
    """Masses [kg], volumes [m3], EuDetCode composition and density, and
    the solver quantities of one trial."""
    mat, wt, d, mb, mmix, p0, T0 = TRIALS[trial]
    mb, mmix = mb*1e-3, mmix*1e-3
    mpart = mmix*wt/100
    mgel = mmix - mpart
    m = {k: mgel*f for k, f in GEL.items()}
    for k, f in C4.items():
        m[k] = mb*f
    V = (sum(m[k]/RHO[k] for k in GEL) + mb/RHO["C4"]
         + (mpart/RHO[mat] if mpart else 0.0))
    # EuDetCode composition (mass fractions without the microballoons)
    comp = {k: v for k, v in m.items() if k != "GMB"}
    if mpart:
        comp["Al"] = mpart
    mfit = sum(comp.values())
    comp = {k: v/mfit for k, v in comp.items()}
    rho_fit = mfit/V
    inert_frac = 0.9 if mat == "Al" else 1.0
    m_inert = mpart*inert_frac
    rho_part = RHO[mat] if mpart else RHO["Al"]
    alpha_part = (m_inert/rho_part)/V if mpart else 0.0
    return dict(trial=trial, mat=mat, wt=wt, d=d, mb=mb, mmix=mmix,
                mtotal=mmix + mb, V=V, R=(3*V/(4*math.pi))**(1/3),
                comp=comp, rho_fit=rho_fit, inert_frac=inert_frac,
                m_inert=m_inert, alpha_part=alpha_part, rho_part=rho_part,
                p0=p0*1e3, T0=T0 + 273.15)


def comp_string(c):
    return " ".join(f"{k}:{v:.6f}" for k, v in c["comp"].items())


def run(trial, eos, outroot):
    c = charge(trial)
    out = os.path.join(outroot, f"trial{trial}")
    exe = os.path.join(EUDET, "build", "src", "CJDetonation")
    with tempfile.TemporaryDirectory() as tmp:
        shutil.copytree(os.path.join(EUDET, "Data"), os.path.join(tmp, "Data"))
        p = os.path.join(tmp, "Data", "Ingredients.json")
        with open(p) as f:
            ing = json.load(f)
        ing["PMMA"] = PMMA_ENTRY
        with open(p, "w") as f:
            json.dump(ing, f, indent=2)
        cmd = [exe, comp_string(c), f"{c['rho_fit']:.2f}", eos,
               "--jwl-products-only", "-o", os.path.abspath(out)]
        if c["mat"] != "baseline":
            cmd[3:3] = ["--inert", f"Al(s):{c['inert_frac']:g}"]
        os.makedirs(out, exist_ok=True)
        with open(os.path.join(out, "log.txt"), "w") as log:
            subprocess.run(cmd, cwd=tmp, stdout=log, stderr=subprocess.STDOUT,
                           check=True)
    with open(os.path.join(out, "command.txt"), "w") as f:
        f.write(" ".join(f'"{a}"' if " " in a else a for a in cmd) + "\n")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--run", nargs="*", default=None, help="trials to fit")
    ap.add_argument("--eos", default="BKW")
    ap.add_argument("--outdir", default=os.path.join(HERE, "eudet_pontalier"))
    a = ap.parse_args()
    print("trial material  wt%  d_um  m_total_kg  R_mm   rho_fit  alpha_part  m_inert_kg")
    for t in TRIALS:
        c = charge(t)
        print(f"{t:>5} {c['mat']:<8} {c['wt']:>4} {c['d']:>5} {c['mtotal']:>10.4f} "
              f"{c['R']*1e3:6.2f} {c['rho_fit']:8.1f} {c['alpha_part']:10.5f} {c['m_inert']:10.4f}")
    for t in a.run or []:
        if t not in TRIALS:
            sys.exit(f"unknown trial {t}")
        print(f"trial {t}: {comp_string(charge(t))} -> {run(t, a.eos, a.outdir)}")


if __name__ == "__main__":
    main()
