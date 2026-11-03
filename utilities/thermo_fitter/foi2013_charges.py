#!/usr/bin/env python3
"""
foi2013_charges.py: charges of the FOI June 2013 TNT and TNT/Al shots
(Fedina, "TNT/Aluminium afterburning in air blasts", FOI-R--3913--SE, 2014,
Table 4.1) and their EuDetCode products-only fits.

Charges: moulded TNT spheres (TNT/Al 80/20 by mass with A100 Al, 42 um),
30 g booster of plastic explosive m/46 (PETN based) in a cavity at the top,
hung at the height of burst. Model conventions (as pontalier_charges.py and
the Comp-B/Al cases):
- the tabulated mass is the main charge; the booster is lumped into the
  composition as PETN/mineral oil 86/14 (m/46, assumed) at 1500 kg/m3;
- TNT at 1630 kg/m3 (as tests/1D_spherical/tnt_blast);
- Al: 90% inert in the fit (--inert "Al(s):0.9"), the separate Al phase.
The ambient state was not recorded: 101325 Pa, 293.15 K.

    python3 foi2013_charges.py [--run SHOT ...] [--eos BKW]
"""

import argparse
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BLAST = os.environ.get("BLAST_DIR", os.path.normpath(os.path.join(HERE, "../..")))
EUDET = os.path.join(BLAST, "external", "eudetcode")

RHO = {"TNT": 1630.0, "booster": 1500.0, "Al": 2700.0}
BOOSTER = {"PETN": 0.86, "Mineral-Oil": 0.14}
P0, T0 = 101325.0, 293.15

# Table 4.1: shot -> (charge, main charge mass [kg], HOB [m])
SHOTS = {
    "1": ("TNT", 0.9103, 0.15), "2": ("TNT", 0.9425, 0.5), "3": ("TNT", 0.9093, 1.0),
    "4": ("TNT", 0.9002, 2.0), "5": ("TNT/Al", 1.0514, 2.0), "6": ("TNT/Al", 1.0146, 1.0),
    "7": ("TNT/Al", 1.0666, 0.5), "8": ("TNT/Al", 1.0009, 0.15), "9": ("TNT", 0.913, 0.5),
    "10": ("TNT", 1.845, 1.0), "11": ("TNT", 1.87, 0.5),
}


def charge(shot):
    """Charge dict in the form used by setup_pontalier.build_case()."""
    kind, m, hob = SHOTS[shot]
    mb = 0.045 if m > 1.5 else 0.030
    wt = 20 if kind == "TNT/Al" else 0
    mAl = m*wt/100
    mTNT = m - mAl
    V = mTNT/RHO["TNT"] + mb/RHO["booster"] + mAl/RHO["Al"]
    comp = {"TNT": mTNT, **{k: mb*f for k, f in BOOSTER.items()}}
    if mAl:
        comp["Al"] = mAl
    mtot = sum(comp.values())
    comp = {k: v/mtot for k, v in comp.items()}
    m_inert = 0.9*mAl
    return dict(shot=shot, kind=kind, hob=hob, mat="Al" if mAl else "baseline", wt=wt,
                d=42 if mAl else 0, mtotal=mtot, V=V, R=(3*V/(4*math.pi))**(1/3),
                comp=comp, rho_fit=mtot/V, m_inert=m_inert,
                alpha_part=(m_inert/RHO["Al"])/V, rho_part=RHO["Al"], p0=P0, T0=T0)


def run(shot, eos, outroot):
    c = charge(shot)
    out = os.path.abspath(os.path.join(outroot, f"shot{shot}"))
    os.makedirs(out, exist_ok=True)
    cmd = [os.path.join(EUDET, "build", "src", "CJDetonation"),
           " ".join(f"{k}:{v:.6f}" for k, v in c["comp"].items()), f"{c['rho_fit']:.2f}", eos]
    if c["mat"] == "Al":
        cmd += ["--inert", "Al(s):0.9"]
    cmd += ["--jwl-products-only", "-o", out]
    with open(os.path.join(out, "log.txt"), "w") as log:
        subprocess.run(cmd, cwd=EUDET, stdout=log, stderr=subprocess.STDOUT, check=True)
    with open(os.path.join(out, "command.txt"), "w") as f:
        f.write(" ".join(f'"{a}"' if " " in a else a for a in cmd) + "\n")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--run", nargs="*", default=None)
    ap.add_argument("--eos", default="BKW")
    ap.add_argument("--outdir", default=os.path.join(HERE, "eudet_foi2013"))
    a = ap.parse_args()
    print("shot charge  m_main  HOB  m_total  R_mm  rho_fit  alpha_Al")
    for s in SHOTS:
        c = charge(s)
        print(f"{s:>4} {c['kind']:<7}{SHOTS[s][1]:7.4f} {c['hob']:4.2f} {c['mtotal']:7.4f} "
              f"{c['R']*1e3:6.2f} {c['rho_fit']:7.1f} {c['alpha_part']:8.5f}")
    for s in a.run or []:
        if s not in SHOTS:
            sys.exit(f"unknown shot {s}")
        print(f"shot {s} -> {run(s, a.eos, a.outdir)}")


if __name__ == "__main__":
    main()
