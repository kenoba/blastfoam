#!/usr/bin/env python3
"""
setup_foi2013_2d.py: 2D axisymmetric case with the ground of one FOI June
2013 shot, started from the 1D free-air stage (../FOI2013_TNT_1D_spherical)
mapped at t0, before the shock reaches the ground.

    python3 setup_foi2013_2d.py --oneD <1D case> --out <case dir>
        [--t0 auto] [--W 7.04] [--H 5.12] [--h 0.032] [--minLevel 2]
        [--mixingLevel 3] [--pressureDelta 2000] [--endTime 18e-3] [--np 6]
    cd <case dir> && ./Allrun

- Mesh: 5-degree wedge about the vertical axis y, x in [0, W], y in [0, H],
  base cells h; the ground y = 0 is a symmetry plane (rigid, inviscid), the
  side and top are transmissive. The charge centre is on the axis at the
  height of burst (from charge.txt of the 1D case).
- Mapping time t0 (auto): T0_BY_HOB, the 1D shock then at about 75-90% of
  the HOB; checked against the 1D probe at r = HOB.
- Physics: constant/ of the 1D case (the shot's products fit, closure,
  transfer and gas chemistry); numerics, AMR and the initial-field set of
  ../ChargeCompBAl_2D_axi_2m (G2b): refinement at 20% jumps of rho, and at
  absolute pressure jumps above --pressureDelta (2 kPa per face, coarsened
  below half of it; 0 keeps the 20% jumps of p of G2). The relative criterion
  does not flag the weak shock beyond about 2.5 m (a 20-30 kPa shock smeared
  over 2-3 cells), which then runs on the base cells (32 mm) and loses 20% of
  its peak (shot 4, 2026-11-11).
  The finest level follows the 1D cell size at the shock radius, as in G2:
  a time table from the 1D probe arrival times (the 1D run must have passed
  its outer probes), never below minLevel (2: 8 mm with 32 mm base cells). The products/air mixing layer is held at
  `mixingLevel` (3: 4 mm), marked by the carbon / (carbon + nitrogen)
  element ratio of the gas, whose products value is computed from the
  transfer composition of the 1D case.
- Allrun decomposes, maps the 1D solution at t0 with rotateFields -refine
  -extend in parallel (centre at the height of burst) and runs blastEulerFoam
  with load balancing; beginTime 0 makes it a restart.
- Probes: ground gauges P1-P3 at 2, 4, 6 m and every 0.5 m along the ground,
  rays up, sideways and down from the charge at the 1D probe radii.
"""

import argparse
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
G2 = os.path.normpath(os.path.join(HERE, "../ChargeCompBAl_2D_axi_2m"))

# Carbon mass fraction of the gas species (as the G2 mixingLayer criterion)
CARBON = {"CO": 0.42882, "CO2": 0.27292, "CH4": 0.74868, "C2H6": 0.79890,
          "C(gr)": 1.0, "CH2O": 0.40002, "CH3": 0.79887, "CH3O": 0.38703,
          "HCO": 0.41392}
NITROGEN = {"N2": 1.0, "N": 1.0, "NH3": 0.82244, "NO": 0.46680}


# Mapping time by height of burst: 1D shock radius about 0.75-0.9 HOB (1D
# shot 4 arrivals: 0.1 m 15 us, 0.3 m 69 us, 0.5 m 154 us, 0.75 m 306 us,
# 1.0 m 515 us)
T0_BY_HOB = {0.15: "1.5e-05", 0.5: "0.0001", 1.0: "0.0004", 2.0: "0.0015"}


def probe_series(case):
    """Probe radii and pressure history of a 1D case (restarts concatenated)."""
    import numpy as np
    base = os.path.join(case, "postProcessing", "pressureProbes")
    dirs = sorted(os.listdir(base), key=float)
    radii, data = None, []
    for d in dirs:
        path = os.path.join(base, d, "p")
        if radii is None:
            radii = [float(l.split("(")[1].split()[0]) for l in open(path) if l.startswith("# Probe")]
        a = np.loadtxt(path, ndmin=2)
        if data:
            a = a[a[:, 0] > data[-1][-1, 0]]
        data.append(a)
    return radii, np.vstack(data)


def arrivals(case, p0):
    """1D arrival time (5% overpressure) at each probe radius, None if not reached."""
    radii, d = probe_series(case)
    out = {}
    for i, r in enumerate(radii):
        hit = d[:, i + 1] - p0 > 0.05*p0
        out[r] = float(d[hit.argmax(), 0]) if hit.any() else None
    return out, float(d[-1, 0])


def cell_size_1d(case):
    """Cell size of the 1D graded mesh as a function of radius."""
    txt = open(os.path.join(case, "system", "blockMeshDict")).read()
    g = re.search(r"simpleGrading\s*\(\(\(([^)]*)\)\s*\(([^)]*)\)\)", txt)
    L1, n1, _ = (float(v) for v in g.group(1).split())
    L2, n2, ratio = (float(v) for v in g.group(2).split())
    q = ratio**(1/(n2 - 1))
    h0 = L2*(q - 1)/(q**n2 - 1)

    def dx(r):
        if r <= L1:
            return L1/n1
        import math
        k = math.log(1 + (r - L1)*(q - 1)/h0)/math.log(q)
        return h0*q**min(k, n2 - 1)
    return dx


def level_table(case, p0, t0, h, minLevel):
    """Finest level (h/2^L about the 1D cell size at the shock radius)
    against time, from the 1D arrivals; a step function."""
    import math
    dx = cell_size_1d(case)
    ta, tend = arrivals(case, p0)
    lev = lambda r: max(minLevel, round(math.log2(h/dx(r))))
    # shock radius at t0: last probe reached by t0
    reached = [r for r, t in ta.items() if t is not None and t <= t0]
    r0 = max(reached) if reached else min(ta)
    later = sorted((t, r) for r, t in ta.items() if t is not None and t > t0)
    reachedAll = [r for r, t in ta.items() if t is not None]
    missing = [r for r, t in ta.items() if t is None]
    drop = [r for r in missing if lev(r) < lev(max(reachedAll))]
    if drop:
        sys.exit(f"1D case at {tend:g} s has not reached r = {min(drop):g} m, where the "
                 f"level drops to {lev(min(drop))}: run it further")
    rows, cur = [(0.0, lev(r0))], lev(r0)
    for t, r in later:
        L = lev(r)
        if L < cur:
            rows += [(t, cur), (t*1.000001, L)]
            cur = L
    rows.append((1.0, cur))
    return rows, r0


def products_ratio(fvModels):
    """C / (C + N) element mass ratio of the transferred products."""
    txt = open(fvModels).read()
    blk = re.search(r"targetFractions\s*\{([^}]*)\}", txt).group(1)
    Y = {k: float(v) for k, v in re.findall(r"(\S+)\s+([0-9.eE+-]+)\s*;", blk)}
    c = sum(Y.get(k, 0)*w for k, w in CARBON.items())
    n = sum(Y.get(k, 0)*w for k, w in NITROGEN.items())
    return c/(c + n)


HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  12
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
FOOTER = "\n// ************************************************************************* //\n"


def block_mesh(W, H, h):
    nx, ny = round(W/h), round(H/h)
    return HEADER.format(obj="blockMeshDict") + f"""\
// Axisymmetric (axis y, vertical) 5-degree wedge, x in [0, {W:g}], y in
// [0, {H:g}] m, base cells {h*1e3:g} mm; ground y = 0 (symmetry plane), side
// and top transmissive.

convertToMeters 1;

W       {W:g};
H       {H:g};
nx      {nx};
ny      {ny};

xw      #calc "$W*cos(2.5*Foam::constant::mathematical::pi/180.0)";
zw      #calc "$W*sin(2.5*Foam::constant::mathematical::pi/180.0)";
nzw     #neg  $zw;

vertices
(
    (0   0  0)
    ($xw 0  $nzw)
    ($xw $H $nzw)
    (0   $H 0)
    ($xw 0  $zw)
    ($xw $H $zw)
);

blocks
(
    hex (0 1 2 3 0 4 5 3) ($nx $ny 1) simpleGrading (1 1 1)
);

boundary
(
    ground
    {{
        type symmetryPlane;
        faces ((0 1 4 0));
    }}
    outlet
    {{
        type patch;
        faces ((1 2 5 4) (3 5 2 3));
    }}
    wedge0
    {{
        type wedge;
        faces ((0 3 2 1));
    }}
    wedge1
    {{
        type wedge;
        faces ((0 4 5 3));
    }}
    // Exposed internal faces during load balancing (redistributor)
    internal
    {{
        type internal;
        faces ();
    }}
);
""" + FOOTER


def control_dict(endTime, hob, W, H, probes_ground, probes_rays):
    pl = "\n".join(f"            ({x:g} {y:g} 0) // {c}" for x, y, c in probes_ground + probes_rays)
    return HEADER.format(obj="controlDict") + f"""\
application     blastEulerFoam;

startFrom       latestTime;

startTime       0;

// The run starts from the mapped 1D solution (Allrun); beginTime 0 makes it
// a restart, so the detonation energy is not added again
beginTime       0;

stopAt          endTime;

endTime         {endTime:g};

writeControl    adjustableRunTime;

writeInterval   2.5e-4;

writeFormat     ascii;

writePrecision  6;

writeCompression off;

timeFormat      general;

timePrecision   6;

runTimeModifiable true;

adjustTimeStep  yes;

deltaT          1e-8;

maxCo           0.2;

functions
{{
    sampleLine
    {{
        type            sets;
        libs            (sampling);
        interpolationScheme cell;
        setFormat       raw;
        executeControl  writeTime;
        writeControl    writeTime;
        sets
        (
            ground
            {{
                type    lineUniform;
                axis    x;
                start   (1e-06 1e-03 0);
                end     ({W - 1e-3:g} 1e-03 0);
                nPoints 3000;
            }}
            axis
            {{
                type    lineUniform;
                axis    y;
                start   (1e-06 1e-06 0);
                end     (1e-06 {H - 1e-3:g} 0);
                nPoints 3000;
            }}
            horizontal
            {{
                type    lineUniform;
                axis    x;
                start   (1e-06 {hob:g} 0);
                end     ({W - 1e-3:g} {hob:g} 0);
                nPoints 3000;
            }}
        );
        fields          (p T.gas U rho alpha.gas alpha.compBAl alpha.Al CO.gas CO2.gas H2.gas O2.gas);
    }}

    pressureProbes
    {{
        type           blastProbes;
        writeControl   timeStep;
        writeInterval  1;
        adjustLocations yes;
        probeLocations
        (
{pl}
        );
        fields          (p);
    }}

    groundPeakPressure
    {{
        type            temporalMax;
        libs            ("libblastFunctionObjects.so");
        field           p;
        writeControl    writeTime;
    }}
}}
""" + FOOTER


def dynamic_mesh(levels, mixingLevel, ratio, pressureDelta):
    src = open(os.path.join(G2, "constant", "dynamicMeshDict")).read()
    body = src[src.index("maxLevel"):]
    if pressureDelta > 0:
        new = ("pressure\n        {\n"
               "            errorEstimator      delta;\n"
               "            field               p;\n"
               f"            lowerRefineLevel    {pressureDelta:g};\n"
               f"            unrefineLevel       {pressureDelta/2:g};\n"
               "            maxRefinement       $maxLevel;\n        }")
        body, n = re.subn(r"pressure\s*\{[^}]*\}", new, body, count=1)
        assert n == 1, "pressure estimator not found in the G2 dynamicMeshDict"
    tab = "\n".join(f"    ({t:<16.10g} {L + 0.5:g})" for t, L in levels)
    body = re.sub(r"maxLevel\s+table\s*\([^;]*?\);",
                  f"maxLevel        table\n(\n{tab}\n);", body, count=1, flags=re.S)
    body = re.sub(r"// Level of the products/air mixing layer[^\n]*\n//[^\n]*\nmixingLevel\s+\d+;",
                  f"// Level of the products/air mixing layer\nmixingLevel     {mixingLevel};", body)
    up = round(ratio - 0.01, 3)
    body = re.sub(r"upperRefineLevel\s+[0-9.]+;", f"upperRefineLevel    {up:g};", body)
    body = re.sub(r"upperUnrefineLevel\s+[0-9.]+;", f"upperUnrefineLevel  {round(ratio - 0.004, 3):g};", body)
    body = re.sub(r"0\.449 in the\s*\n(\s*)// products", f"{ratio:.3f} in the\n\\1// products", body)
    body = body.replace("Refined between 0.01 and 0.44", f"Refined between 0.01 and {up:g}")
    head = (HEADER.format(obj="dynamicMeshDict")
            + "// FOI 2013 2D ground case (setup_foi2013_2d.py): AMR of ChargeCompBAl_2D_axi_2m\n"
            + "// (G2b). Finest level (floored) against time: about the 1D cell size at the\n"
            + "// shock radius (from the 1D arrival times), at least minLevel;\n"
            + f"// mixing layer at level {mixingLevel}; products C/(C+N) element ratio\n"
            + (f"// Pressure: absolute jumps above {pressureDelta:g} Pa per face (weak shocks).\n"
               if pressureDelta > 0 else "")
            + f"// {ratio:.4f} (1D transfer composition).\n\n")
    return head + body


ALLRUN = """#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Mapping time and the 1D free-air case (oneD -> ../FOI2013_TNT_1D_spherical run)
t0=__T0__
centre='(0 __HOB__ 0)'
oneD=oneD
if [ ! -d "$oneD/$t0" ]
then
    echo "no 1D solution at $oneD/$t0" >&2
    exit 1
fi

rm -rf [1-9]* 0.[0-9]* postProcessing log.* constant/polyMesh dynamicCode processor*

runApplication blockMesh
mkdir "$t0"
for f in 0/*.orig
do
    cp "$f" "$t0/$(basename "$f" .orig)"
done

# Map in parallel after decomposePar (refinement history, see
# ../ChargeCompBAl_2D_axi_2m); all 1D fields except derived and
# post-processing ones
runApplication decomposePar
fields=$(ls "$oneD/$t0" | grep -v -x -e uniform -e pTemporalMax -e p.compBAl \\
    -e p.gas -e p.Al -e Ptot.Al -e Qdot.gas -e V -e Vc -e C -e Ccx -e Ccy \\
    -e Ccz | tr '\\n' ' ')
runParallel rotateFields "$(cd "$oneD" && pwd -P)" -sourceTime "$t0" \\
    -sourceCentre '(0 0 0)' -centre "$centre" \\
    -additionalFields "($fields)" -refine -extend
if ! grep -q "^Finished" log.rotateFields
then
    echo "rotateFields failed: see log.rotateFields" >&2
    exit 1
fi

touch case.foam
runParallel $(getApplication)
"""

ALLCLEAN = """#!/bin/sh
cd ${0%/*} || exit 1    # run from this directory
rm -rf [1-9]* 0.[0-9]* postProcessing log.* constant/polyMesh dynamicCode processor* *.foam
for f in 0/*; do
    case "$f" in
        *.orig) ;;
        *) rm -rf "${f:?}" ;;
    esac
done
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--oneD", required=True, help="1D case (setup_foi2013.py) run past t0")
    ap.add_argument("--out", required=True)
    ap.add_argument("--t0", default="auto", help="mapping time, as the 1D time directory name")
    ap.add_argument("--W", type=float, default=7.04)
    ap.add_argument("--H", type=float, default=5.12)
    ap.add_argument("--h", type=float, default=0.032)
    ap.add_argument("--minLevel", type=int, default=2)
    ap.add_argument("--mixingLevel", type=int, default=3)
    ap.add_argument("--endTime", type=float, default=18e-3)
    ap.add_argument("--pressureDelta", type=float, default=2000.0,
                    help="refine at absolute pressure jumps above this [Pa]; 0: 20%% jumps")
    ap.add_argument("--np", type=int, default=6)
    a = ap.parse_args()

    oneD = os.path.abspath(a.oneD)
    info = open(os.path.join(oneD, "charge.txt")).read().split()
    kv = dict(zip(info[0::2], info[1::2]))
    hob = float(kv["hob"])
    p0 = float(kv["p0"])
    t0 = T0_BY_HOB[hob] if a.t0 == "auto" else a.t0
    if not os.path.isdir(os.path.join(oneD, t0)):
        sys.exit(f"no 1D solution at {oneD}/{t0}")
    ta, _ = arrivals(oneD, p0)
    rh = min((r for r in ta if r >= hob - 1e-9), default=None)
    if rh is not None and ta[rh] is not None and ta[rh] <= float(t0):
        sys.exit(f"the 1D shock reaches r = {rh:g} m (HOB) at {ta[rh]:g} s, before t0 = {t0}")
    levels, r0 = level_table(oneD, p0, float(t0), a.h, a.minLevel)
    if os.path.exists(a.out):
        sys.exit(f"{a.out} exists")
    out = os.path.abspath(a.out)
    os.makedirs(os.path.join(out, "system"))

    # Physics from the 1D case, numerics and fields from G2
    shutil.copytree(os.path.join(oneD, "constant"), os.path.join(out, "constant"),
                    ignore=shutil.ignore_patterns("polyMesh", "dynamicMeshDict"))
    for f in ("fvSchemes", "fvSolution", "decomposeParDict"):
        shutil.copy(os.path.join(G2, "system", f), os.path.join(out, "system", f))
    dp = os.path.join(out, "system", "decomposeParDict")
    s = open(dp).read()
    open(dp, "w").write(re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains  {a.np};", s))
    os.makedirs(os.path.join(out, "0"))
    for f in os.listdir(os.path.join(G2, "0")):
        if f.endswith(".orig"):
            s = open(os.path.join(G2, "0", f)).read()
            s = re.sub(r"^(\s*)symmetry(\s*)$", r"\1ground\2", s, flags=re.M)
            open(os.path.join(out, "0", f), "w").write(s)

    fvm = os.path.join(out, "constant", "fvModels")
    if os.path.isfile(fvm) and "turbulentMixingDiffusion" in open(fvm).read():
        # Tight species solves: the gas thermo normalises the mass fractions after the
        # species are solved one by one, so the diffusion solver error would turn
        # into element changes (turbulentMixingDiffusion warns above 1e-12)
        fs = os.path.join(out, "system/fvSolution")
        s, n = re.subn(r'("[^"\n]*\\.gas[^"\n]*"\s*\{[^}]*?tolerance\s+)[^;]+;', r"\g<1>1e-15;",
                       open(fs).read())
        if n != 1:
            sys.exit(f"{fs}: no single species solver entry to tighten")
        open(fs, "w").write(s)

    ratio = products_ratio(os.path.join(oneD, "constant", "fvModels"))
    open(os.path.join(out, "constant", "dynamicMeshDict"), "w").write(
        dynamic_mesh(levels, a.mixingLevel, ratio, a.pressureDelta))
    open(os.path.join(out, "system", "blockMeshDict"), "w").write(block_mesh(a.W, a.H, a.h))

    ground = [(x, 1e-3, f"ground x {x:g}" + (f" P{int(x/2)}" if x in (2, 4, 6) else ""))
              for x in [0.0 + 1e-6] + [0.5*i for i in range(1, int(a.W/0.5) + 1) if 0.5*i < a.W]]
    rays = []
    for r in (0.5, 1.0, 1.5, 1.9):
        rays += [(1e-6, hob + r, f"up r {r:g}"), (r, hob, f"side r {r:g}")]
        if r < hob:
            rays.append((1e-6, hob - r, f"down r {r:g}"))
    for r in (2.5, 3.0):
        if hob + r < a.H:
            rays.append((1e-6, hob + r, f"up r {r:g}"))
        rays.append((r, hob, f"side r {r:g}"))
    open(os.path.join(out, "system", "controlDict"), "w").write(
        control_dict(a.endTime, hob, a.W, a.H, ground, rays))

    for name, txt in (("Allrun", ALLRUN.replace("__T0__", t0).replace("__HOB__", f"{hob:g}")),
                      ("Allclean", ALLCLEAN)):
        p = os.path.join(out, name)
        open(p, "w").write(txt)
        os.chmod(p, 0o755)
    os.symlink(oneD, os.path.join(out, "oneD"))
    with open(os.path.join(out, "charge.txt"), "w") as f:
        f.write(" ".join(info) + f" t0 {t0} W {a.W:g} H {a.H:g} h {a.h:g} minLevel {a.minLevel}"
                f" mixingLevel {a.mixingLevel}"
                f" pressureDelta {a.pressureDelta:g}\n")
    print(f"{out}: HOB {hob:g} m, products C/(C+N) {ratio:.4f}, "
          f"{round(a.W/a.h)} x {round(a.H/a.h)} base cells, t0 {t0} (1D shock past r = {r0:g} m), "
          f"levels {[(round(t, 6), L) for t, L in levels]}")


if __name__ == "__main__":
    main()
