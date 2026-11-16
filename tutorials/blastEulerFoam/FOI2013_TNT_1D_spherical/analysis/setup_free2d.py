#!/usr/bin/env python3
"""2D axisymmetric free-field FOI shot 4 run from the detonation (no ground,
no 1D mapping), with the initiation point on the axis at height yi above
the charge centre (booster at the top: yi ~ R - r_booster).
    setup_free2d.py <1D case> <out> <yi> [--chem off] [--np 6] [--end 6.5e-3]"""
import argparse, os, re, shutil, sys, math
T2D = "/opt/OpenFOAM/blastfoam/tutorials/blastEulerFoam/FOI2013_TNT_2D_axi_ground"
sys.path.insert(0, T2D)
import setup_foi2013_2d as s
ap = argparse.ArgumentParser()
ap.add_argument("oneD"); ap.add_argument("out"); ap.add_argument("yi", type=float)
ap.add_argument("--chem", default="off"); ap.add_argument("--np", type=int, default=6)
ap.add_argument("--end", type=float, default=6.5e-3); ap.add_argument("--L", type=float, default=3.2)
ap.add_argument("--h", type=float, default=0.032); ap.add_argument("--maxLevel", type=int, default=5)
a = ap.parse_args()
oneD = os.path.abspath(a.oneD); out = os.path.abspath(a.out)
if os.path.exists(out): sys.exit(f"{out} exists")
info = open(os.path.join(oneD, "charge.txt")).read().split()
kv = dict(zip(info[0::2], info[1::2])); R = float(kv["R"]); p0 = float(kv["p0"])
os.makedirs(os.path.join(out, "system")); os.makedirs(os.path.join(out, "0"))
shutil.copytree(os.path.join(oneD, "constant"), os.path.join(out, "constant"),
                ignore=shutil.ignore_patterns("polyMesh", "dynamicMeshDict"))
# initiation point
f = os.path.join(out, "constant", "physicalProperties.compBAl"); t = open(f).read()
t, n = re.subn(r"points\s*\(\s*\(0 0 0\)\s*\);", f"points ((0 {a.yi:g} 0));", t); assert n == 1
open(f, "w").write(t)
if a.chem == "off":
    for c in ("chemistryProperties", "chemistryProperties.gas"):
        f = os.path.join(out, "constant", c); t = open(f).read()
        open(f, "w").write(re.sub(r"^chemistry\s+on;", "chemistry           off;", t, flags=re.M))
for f in ("fvSchemes", "fvSolution", "decomposeParDict"):
    shutil.copy(os.path.join(s.G2, "system", f), os.path.join(out, "system", f))
dp = os.path.join(out, "system", "decomposeParDict")
txt = open(dp).read()
open(dp, "w").write(re.sub(r"numberOfSubdomains\s+\d+;", f"numberOfSubdomains  {a.np};", txt))
for f in os.listdir(os.path.join(s.G2, "0")):
    if f.endswith(".orig"):
        t = open(os.path.join(s.G2, "0", f)).read()
        t = re.sub(r"^(\s*)symmetry(\s*)$", r"\1bottom\2", t, flags=re.M)
        open(os.path.join(out, "0", f), "w").write(t)
# levels from the 1D cell size at the shock, capped at maxLevel
lv, _ = s.level_table(oneD, p0, 0.0, a.h, 2)
lv = [(t_, min(L, a.maxLevel)) for t_, L in lv]
rows = [lv[0]]
for t_, L in lv[1:]:
    if L != rows[-1][1]:
        rows += [(t_*(1 - 2e-5), rows[-1][1]), (t_, L)]
lv = rows + [(1.0, rows[-1][1])]
ratio = s.products_ratio(os.path.join(oneD, "constant", "fvModels"))
open(os.path.join(out, "constant", "dynamicMeshDict"), "w").write(s.dynamic_mesh(lv, 3, ratio, 2000.0))
# mesh: wedge about y, x in [0, L], y in [-L, L]; all outer faces transmissive
bm = s.block_mesh(a.L, 2*a.L, a.h).replace(f"H       {2*a.L:g};", f"H       {a.L:g};")
bm = bm.replace("(0   0  0)", f"(0   {-a.L:g}  0)").replace("($xw 0  $nzw)", f"($xw {-a.L:g}  $nzw)") \
       .replace("($xw 0  $zw)", f"($xw {-a.L:g}  $zw)")
bm = re.sub(r"ground\s*\{\s*type symmetryPlane;\s*faces \(\(0 1 4 0\)\);\s*\}",
            "bottom\n    {\n        type patch;\n        faces ((0 1 4 0));\n    }", bm)
open(os.path.join(out, "system", "blockMeshDict"), "w").write(bm)
# 0/*.orig: the 'bottom' patch as the outlet
for f in os.listdir(os.path.join(out, "0")):
    p = os.path.join(out, "0", f); t = open(p).read()
    m = re.search(r"(\n\s*outlet\s*\{[^}]*\})", t)
    t = re.sub(r"\n\s*bottom\s*\{[^}]*\}", "", t)
    if m: t = t.replace(m.group(1), m.group(1) + m.group(1).replace("outlet", "bottom", 1))
    open(p, "w").write(t)
# probes: rays at elevation angles (deg) from the charge centre
probes = []
for ang in (90, 0, -18.43, -26.57, -45, -90):
    for r in (0.5, 1.0, 1.5, 2.0, 2.5, 2.828):
        x, y = r*math.cos(math.radians(ang)), r*math.sin(math.radians(ang))
        probes.append((max(x, 1e-6), y, f"ang {ang:g} r {r:g}"))
cd = s.control_dict(a.end, 0.0, a.L, a.L, probes, [])
cd = cd.replace("beginTime       0;", "").replace("// The run starts from the mapped 1D solution (Allrun); beginTime 0 makes it\n// a restart, so the detonation energy is not added again\n", "")
cd = re.sub(r"ground\s*\{[^}]*\}", "", cd, count=1)
cd = re.sub(r"axis\s*\{\s*type    lineUniform;\s*axis    y;\s*start   \(1e-06 1e-06 0\);",
            f"axis\n            {{\n                type    lineUniform;\n                axis    y;\n                start   (1e-06 {-a.L + 1e-3:g} 0);", cd)
cd = cd.replace("groundPeakPressure", "peakPressure")
open(os.path.join(out, "system", "controlDict"), "w").write(cd)
open(os.path.join(out, "system", "setFieldsDict"), "w").write(s.HEADER.format(obj="setFieldsDict") + f"""
fields (alpha.compBAl alpha.Al alpha.gas);
nBufferLayers 2;
defaultFieldValues
(
    volScalarFieldValue alpha.compBAl 0
    volScalarFieldValue alpha.Al 0
    volScalarFieldValue alpha.gas 1
    volScalarFieldValue lambda.compBAl 0
);
regions
(
    sphereToCell
    {{
        centre (0 0 0);
        radius {R:g};
        level {a.maxLevel};
        refineInternal yes;
        backup {{ centre (0 0 0); radius {R*1.2:g}; }}
        fieldValues
        (
            volScalarFieldValue alpha.compBAl 1
            volScalarFieldValue alpha.Al 0
            volScalarFieldValue alpha.gas 0
            volScalarFieldValue lambda.compBAl 0
        );
    }}
);
""" + s.FOOTER)
allrun = """#!/bin/sh
cd ${0%/*} || exit 1
. $WM_PROJECT_DIR/bin/tools/RunFunctions
for f in 0/*.orig; do cp "$f" "0/$(basename "$f" .orig)"; done
runApplication blockMesh
# refine after decomposePar (refinement history)
runApplication decomposePar
runParallel setRefinedFields
touch case.foam
runParallel $(getApplication)
"""
open(os.path.join(out, "Allrun"), "w").write(allrun); os.chmod(os.path.join(out, "Allrun"), 0o755)
open(os.path.join(out, "charge.txt"), "w").write(" ".join(info) + f" yi {a.yi:g} chem {a.chem} L {a.L:g} h {a.h:g} maxLevel {a.maxLevel}\n")
print(out, "levels", [(round(t_, 6), L) for t_, L in lv], "R", R)
