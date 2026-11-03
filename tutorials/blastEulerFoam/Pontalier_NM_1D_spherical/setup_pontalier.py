#!/usr/bin/env python3
"""
setup_pontalier.py: 1D spherical case of one Pontalier et al. (2021) trial
(gelled nitromethane with Al or glass particles, Shock Waves 31:203-230).

    python3 setup_pontalier.py --trial 6 --out <case dir>
        [--closure laminar|pasrd] [--R 4.0] [--endTime 4e-3]

The case is built from tutorials/blastEulerFoam/ChargeCompBAl_1D_spherical_2m
with setup_e2.py (Phase D Al models, products -> gas transfer, gas chemistry):
variant full for Al and the baseline, products (Al never ignites) for glass.
The explosive phase keeps the name compBAl. Per trial:

- products: EuDetCode products-only fit of
  utilities/thermo_fitter/eudet_pontalier/trial<N> (pontalier_charges.py;
  booster lumped into the charge), entries from products_only_fit.py on the
  gas thermo of this case; products Cv from the transfer gas at 2500 K;
  reactants Murnaghan with the EuDetCode NM data (K0 5 GPa, n 6);
- charge radius, particle volume fraction and density from
  pontalier_charges.py; particle diameter from Table 1; glass beads are the
  Al phase with glass properties (2500 kg/m3, Cv 840, no ignition);
- ambient p0 and T0 of the trial (altitude 700 m);
- domain to R (default 4 m: 200 cells of 0.5 mm to 0.1 m, then graded to
  15 mm), probes at the gauges (1.0, 1.5, 2.0 m) and every 0.5 m to 3.5 m;
- closure pasrd: PaSR with the cell-size mixing time (G3c).
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
BLAST = os.environ.get("BLAST_DIR", os.path.normpath(os.path.join(HERE, "../../..")))
BASE = os.path.join(BLAST, "tutorials/blastEulerFoam/ChargeCompBAl_1D_spherical_2m")
FIT = os.path.join(BLAST, "utilities/thermo_fitter")
sys.path.insert(0, FIT)
import pontalier_charges as pc  # noqa: E402
import products_only_fit as pof  # noqa: E402

R_AIR = 287.1  # J/(kg K): 1.225 kg/m3 at 101298 Pa, 288 K in the base case
CV_REACTANTS = 1030.0


def sub(path, pattern, repl, count=1, flags=0):
    with open(path) as f:
        s = f.read()
    s, n = re.subn(pattern, repl, s, count=count, flags=flags)
    if n != count:
        sys.exit(f"{path}: {n} substitutions of {pattern!r}, expected {count}")
    with open(path, "w") as f:
        f.write(s)


def set_uniform(case, field, value):
    sub(os.path.join(case, "0", field + ".orig"),
        r"internalField\s+uniform\s+[^;]+;", f"internalField   uniform {value};")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--trial", required=True, choices=list(pc.TRIALS))
    ap.add_argument("--out", required=True)
    ap.add_argument("--closure", choices=("laminar", "pasrd"), default="laminar")
    ap.add_argument("--R", type=float, default=4.0)
    ap.add_argument("--endTime", type=float, default=4e-3)
    a = ap.parse_args()

    eudet = os.path.join(FIT, "eudet_pontalier", f"trial{a.trial}")
    if not os.path.isfile(os.path.join(eudet, "cj_output.txt")):
        sys.exit(f"no fit in {eudet}: run pontalier_charges.py --run {a.trial}")
    c = pc.charge(a.trial)
    label = (f"Pontalier et al. (2021) trial {a.trial} ({c['mat']}, {c['wt']}%): "
             f"EuDetCode products-only fit eudet_pontalier/trial{a.trial}")
    build_case(a.out, eudet, c, label, closure=a.closure, R=a.R, endTime=a.endTime,
               info=f"trial {a.trial}")


def build_case(out, eudet, c, label, closure="laminar", R=4.0, endTime=4e-3,
               K0=5e9, n_murnaghan=6.0, probes=None, info=""):
    """Case from the 2 m Comp-B/Al tutorial for an explosive with the
    EuDetCode products-only fit in eudet and the charge c (keys mat, wt, d,
    R, alpha_part, p0, T0, mtotal; mat 'baseline' for no particles)."""
    out = os.path.abspath(out)
    if os.path.exists(out):
        sys.exit(f"{out} exists")

    shutil.copytree(BASE, out, ignore=shutil.ignore_patterns(
        "__pycache__", "postProcessing", "log.*", "[1-9]*", "0.[0-9]*",
        "polyMesh", "dynamicCode", "*.foam", "compare_e2.py", "README.md"))
    variant = "products" if c["mat"] == "glass" else "full"
    subprocess.run([sys.executable, os.path.join(out, "setup_e2.py"), "--variant", variant],
                   cwd=out, check=True)
    os.remove(os.path.join(out, "setup_e2.py"))

    # products fit on the gas thermo of this case
    r = pof.compute(eudet, out, T0=c["T0"], cv_reactants=CV_REACTANTS,
                    inert=("Al(s)", "Al"))
    D_ray, D_eos, _ = pof.burn_consistency(r, CV_REACTANTS*c["T0"])
    vDet = r["D"] if r["D"] > D_eos*1.0005 else D_eos*1.0005
    W = r["cj"].get("Molar Mass (W_eq, kg/mol)", r["W_mix"]*1e-3)*1e3

    pp = os.path.join(out, "constant", "physicalProperties.compBAl")
    sub(pp, r"rho0 1624\.45;", f"rho0 {r['rho0']:.6g};", count=2)
    sub(pp, r"n 6\.0;\n    K0 12\.61e9;", f"n {n_murnaghan:g};\n    K0 {K0:g};")
    for k, v in (("A", r["A"]), ("B", r["B"]), ("C", r["C"])):
        sub(pp, rf"\n    {k} [^;]+;", f"\n    {k} {v:.6g};")
    sub(pp, r"R1 [^;]+;", f"R1 {r['R1']:.6g};")
    sub(pp, r"R2 [^;]+;", f"R2 {r['R2']:.6g};")
    sub(pp, r"omega [^;]+;", f"omega {r['w']:.6g};")
    sub(pp, r"Cv 1450;", f"Cv {r['Cv2500']:.1f};", count=2)
    sub(pp, r"T_CJ [^;]+;", f"T_CJ {r['T_CJ']:.6g};")
    sub(pp, r"rho_CJ [^;]+;", f"rho_CJ {r['rho_CJ']:.6g};")
    sub(pp, r"omegaT [^;]+;", f"omegaT {r['omegaT']:.6f};")
    sub(pp, r"E0 [^;]+;", f"E0 {r['E0']:.2f};")
    sub(pp, r"specie \{ molWeight 24\.02; \}", f"specie {{ molWeight {W:.3f}; }}")
    sub(pp, r"//    = [^\n]*J/kg", f"//    = {r['es_CJ']:.5g} - {CV_REACTANTS*c['T0']:.5g} - {r['hug']:.5g} J/kg")
    sub(pp, r"e0 [0-9.e+-]+;", f"e0 {r['e0']:.6g};")
    sub(pp, r"vDet [^;]+;", f"vDet {vDet:.2f};")
    sub(pp, r"Pcj [^;]+;", f"Pcj {r['p_CJ']:.6g};")
    sub(pp, r"type detonating;",
        f"// {label}\n"
        f"// (products_only_fit.py; D {r['D']:.2f} m/s, "
        f"CJ speed of the products EOS {D_eos:.1f} m/s)\ntype detonating;")

    # transfer composition: freeze-out gas of this fit
    fv = os.path.join(out, "constant", "fvModels")
    tf = "".join(f"        {k:<12}{y:.8f};\n" for k, y in r["Y"].items())
    sub(fv, r"    targetFractions\n    \{\n.*?\n    \}", "    targetFractions\n    {\n" + tf + "    }",
        flags=re.S)

    # particles
    if c["mat"] == "glass":
        pa = os.path.join(out, "constant", "physicalProperties.Al")
        sub(pa, r"rho0\s+2700;", "rho0            2500;")
        sub(pa, r"K0\s+72\.7e9;\n    n\s+4\.3;", "K0              37e9;\n    n               4;")
        sub(pa, r"thermodynamics\n\{.*?\n\}", "thermodynamics\n{\n    // soda-lime glass beads\n"
            "    type            eConst;\n    Cv              840;\n    Hf              0;\n}", flags=re.S)
        sub(pa, r"thermo\s+eMelting;", "thermo          eConst;")
        sub(pa, r"kappa\s+237;", "kappa           1;")
        sub(pa, r"residualRho\s+2700;", "residualRho     2500;")
        sub(os.path.join(out, "constant", "phaseProperties"), r"residualRho\s+2700;", "residualRho     2500;")
        set_uniform(out, "rho.Al", 2500)
    set_uniform(out, "d.Al", f"{max(c['d'], 1)*1e-6:g}")

    # charge
    sf = os.path.join(out, "system", "setFieldsDict")
    sub(sf, r"radius 0\.05;", f"radius {c['R']:.5f};")
    sub(sf, r"alpha\.compBAl 0\.818", f"alpha.compBAl {1 - c['alpha_part']:.6f}")
    sub(sf, r"alpha\.Al 0\.182", f"alpha.Al {c['alpha_part']:.6f}")
    set_uniform(out, "rho.compBAl", f"{r['rho0']:.6g}")

    # ambient
    for f in ("T", "T.gas", "T.compBAl", "T.Al"):
        set_uniform(out, f, f"{c['T0']:.2f}")
    set_uniform(out, "p", f"{c['p0']:.0f}")
    set_uniform(out, "rho.gas", f"{c['p0']/(R_AIR*c['T0']):.6g}")

    # mesh: 0.5 mm to 0.1 m, graded to 15 mm at R
    L = R - 0.1
    ratio = 30.0
    n = round(L*math.log(ratio)/(0.015 - 0.0005))
    bm = os.path.join(out, "system", "blockMeshDict")
    sub(bm, r"R       2\.0;", f"R       {R:g};")
    sub(bm, r"\(799 1 1\) simpleGrading \(\(\(0\.1 200 1\) \(1\.9 599 20\)\) 1 1\)",
        f"({200 + n} 1 1) simpleGrading (((0.1 200 1) ({L:g} {n} {ratio:g})) 1 1)")

    # run control and probes
    cd = os.path.join(out, "system", "controlDict")
    sub(cd, r"endTime         2e-3;", f"endTime         {endTime:g};")
    sub(cd, r"writeInterval   2e-5;", "writeInterval   5e-5;")
    sub(cd, r"end     \(2\.0 0 0\);\n                nPoints 2000;",
        f"end     ({R:g} 0 0);\n                nPoints {int(R*1000)};")
    probes = probes or [0.05, 0.1, 0.15, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5]
    sub(cd, r"probeLocations\n        \(\n.*?\n        \);",
        "probeLocations\n        (\n" + "".join(f"            ({x:g} 0 0)\n" for x in probes if x < R)
        + "        );", flags=re.S)

    if closure == "pasrd":
        cp = os.path.join(out, "constant", "combustionProperties.gas")
        sub(cp, r"combustionModel\s+laminar;",
            "combustionModel  PaSR;\n\nPaSRCoeffs\n{\n    Cmix            1;\n"
            "    mixingTime      delta;\n    Cs              0.17;\n}")

    shutil.copy(os.path.join(HERE, "summary.py"), out)
    with open(os.path.join(out, "charge.txt"), "w") as f:
        f.write(f"{info} closure {closure} material {c['mat']} wt {c['wt']} "
                f"d_um {c['d']} m_total {c['mtotal']:.4f} R {c['R']:.5f} "
                f"alpha_part {c['alpha_part']:.6f} p0 {c['p0']:.0f} T0 {c['T0']:.2f} "
                f"vDet {vDet:.2f} pCJ {r['p_CJ']:.4g} rho0 {r['rho0']:.6g} cells {200 + n}\n")
    print(open(os.path.join(out, "charge.txt")).read().strip())
    return r


if __name__ == "__main__":
    main()
