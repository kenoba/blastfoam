#!/usr/bin/env python3
"""
setup_foi2013.py: 1D spherical case of one FOI June 2013 shot (1 kg TNT or
TNT/Al 80/20, Fedina, FOI-R--3913--SE, 2014), the free-air stage before the
shock reaches the ground; its fields are mapped into the 2D axisymmetric
case with the ground.

    python3 setup_foi2013.py --shot 4 --out <case dir>
        [--closure laminar|pasrd] [--R 4.0] [--endTime 4e-3]

Built with build_case() of ../Pontalier_NM_1D_spherical/setup_pontalier.py
(the 2 m Comp-B/Al tutorial with the Phase D Al models, products -> gas
transfer and gas chemistry); the products are the EuDetCode products-only
fit utilities/thermo_fitter/eudet_foi2013/shot<N> (foi2013_charges.py),
the reactants Murnaghan with the EuDetCode TNT data (K0 8.93 GPa, n 5.5).
Extra probes at the slant distances of the ground gauges P1-P3 (2, 4, 6 m
from the point below the charge).
"""

import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BLAST = os.environ.get("BLAST_DIR", os.path.normpath(os.path.join(HERE, "../../..")))
sys.path.insert(0, os.path.join(BLAST, "tutorials/blastEulerFoam/Pontalier_NM_1D_spherical"))
sys.path.insert(0, os.path.join(BLAST, "utilities/thermo_fitter"))
import foi2013_charges as fc  # noqa: E402
from setup_pontalier import build_case  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--shot", required=True, choices=list(fc.SHOTS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--closure", choices=("laminar", "pasrd"), default="laminar")
    ap.add_argument("--R", type=float, default=4.0)
    ap.add_argument("--endTime", type=float, default=4e-3)
    a = ap.parse_args()
    eudet = os.path.join(BLAST, "utilities/thermo_fitter/eudet_foi2013", f"shot{a.shot}")
    if not os.path.isfile(os.path.join(eudet, "cj_output.txt")):
        sys.exit(f"no fit in {eudet}: run foi2013_charges.py --run {a.shot}")
    c = fc.charge(a.shot)
    slant = [round(math.hypot(x, c["hob"]), 3) for x in (2.0, 4.0, 6.0)]
    probes = sorted(set([0.1, 0.15, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5] + slant))
    label = (f"FOI 2013 shot {a.shot} ({c['kind']}, HOB {c['hob']:g} m): "
             f"EuDetCode products-only fit eudet_foi2013/shot{a.shot}")
    build_case(a.out, eudet, c, label, closure=a.closure, R=a.R, endTime=a.endTime,
               K0=8.93e9, n_murnaghan=5.5, probes=probes,
               info=f"shot {a.shot} hob {c['hob']:g}")
    print(f"ground gauges P1-P3 at slant distances {slant} m"
          + ("" if all(s < a.R for s in slant) else f" (beyond R = {a.R:g} m: not probed)"))


if __name__ == "__main__":
    main()
