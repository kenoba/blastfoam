#!/usr/bin/env python3
"""Set up the 1D turbulentMixingDiffusion test (blastEulerFoam, one gas phase).

A closed planar tube, x in [-L, L], at rest and 1 atm. Pure N2 (the "air",
YAir = 1) for x < 0 and pure CO (the "products", YProducts = 0) for x > 0.
N2 and CO have nearly equal molar masses, so mixing at constant temperature
leaves the pressure uniform and the gas at rest, and the mass fraction obeys
dY/dt = D(t) d2Y/dx2 with the solution

    Y_CO = 0.5 erfc(-x/(2 sqrt(S))),  S = int_0^t D dt'.

D(t) = k t (a table Function1), so S = k t^2/2.

usage: setup_case.py [--TLeft 298.15] [--TRight 298.15] [--energy yes|no]
                     [--k 1250] [--endTime 1e-3] [--nCells 200] [--D0 0]
Writes system/, constant/phaseProperties, constant/fvModels and 0/ in the
case directory (the directory of this script, or --case)."""
import argparse
import os

ap = argparse.ArgumentParser()
ap.add_argument("--case", default=os.path.dirname(os.path.abspath(__file__)))
ap.add_argument("--TLeft", type=float, default=298.15)
ap.add_argument("--TRight", type=float, default=298.15)
ap.add_argument("--energy", default="no")
ap.add_argument("--k", type=float, default=1250.0)
ap.add_argument("--endTime", type=float, default=1e-3)
ap.add_argument("--deltaT", type=float, default=1e-6)
ap.add_argument("--nCells", type=int, default=200)
ap.add_argument("--L", type=float, default=0.25)
a = ap.parse_args()

P0 = 101325.0
XO2 = 0.8*a.L
SPECIES = ["O2", "N2", "CO"]

BANNER = """/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  12
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""
END = "\n// ************************************************************************* //\n"


def write(rel, body, cls="dictionary"):
    path = os.path.join(a.case, rel)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    loc, obj = os.path.split(rel)
    with open(path, "w") as f:
        f.write(BANNER.format(cls=cls, loc=loc, obj=obj) + body + END)


N = a.nCells
dx = 2*a.L/N
xc = [-a.L + (i + 0.5)*dx for i in range(N)]

# ---------------------------------------------------------------- system
write("system/blockMeshDict", f"""convertToMeters 1;

vertices
(
    ({-a.L} 0 0) ({a.L} 0 0) ({a.L} 0.01 0) ({-a.L} 0.01 0)
    ({-a.L} 0 0.01) ({a.L} 0 0.01) ({a.L} 0.01 0.01) ({-a.L} 0.01 0.01)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({N} 1 1) simpleGrading (1 1 1)
);

boundary
(
    ends
    {{
        type symmetry;
        faces ((0 4 7 3) (1 2 6 5));
    }}
    sides
    {{
        type empty;
        faces ((0 1 5 4) (3 2 6 7) (0 1 2 3) (4 5 6 7));
    }}
);
""")

write("system/controlDict", f"""application     blastEulerFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {a.endTime:g};

libs            ("libblastThermodynamics.so");

deltaT          {a.deltaT:g};
adjustTimeStep  no;
maxCo           0.5;

writeControl    runTime;
writeInterval   {a.endTime/4:g};
writeFormat     ascii;
writePrecision  12;
writeCompression off;
timeFormat      general;
timePrecision   8;
runTimeModifiable false;
""")

write("system/fvSchemes", """fluxSchemes
{
    gas
    {
        fluxScheme  HLLC;
    }
}

ddtSchemes
{
    default         Euler;
    timeIntegrator  RK2SSP;
}

gradSchemes
{
    default         cellMDLimited leastSquares 1.0;
}

divSchemes
{
    default         none;
    div(alphaRhoPhi.gas,Yi) Riemann;
    div((((alpha.gas*rho.gas)*nuEff.gas)*dev2(T(grad(U.gas))))) Gauss linear;
}

laplacianSchemes
{
    default         Gauss linear corrected;
}

interpolationSchemes
{
    default             linear;
    "reconstruct(.*)"   vanLeer;
}

snGradSchemes
{
    default         corrected;
}
""")

write("system/fvSolution", """solvers
{
    "(rho|alphaRhoU|alphaRhoE|alpha).*"
    {
        solver          diagonal;
    }

    "(U|e|Yi|.*\\.gas|Ydefault).*"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-12;
        relTol          0;
    }
}
""")

# -------------------------------------------------------------- constant
write("constant/phaseProperties", """phases          (gas);

gas
{
    phaseType       fluid;
    diameterModel   constant;
    constantDiameterCoeffs
    {
        d               0;
    }
    residualRho     1e-06;
    residualAlpha   1e-06;
}

blending
{
    default
    {
        type            none;
        continuousPhase gas;
    }
}

interfacialPressure ();
interfacialVelocity ();
massTransfer ();
aspectRatio ();
solveDragODE off;
drag ();
virtualMass ();
heatTransfer ();
lift ();
wallLubrication ();
turbulentDispersion ();

pressureRelaxationSolver
{
    type            none;
}

pressureRelaxation ();
""")

ts = [a.endTime*i/20 for i in range(21)] + [10*a.endTime]
table = "\n".join(f"            ({t:.8e} {a.k*t:.8e})" for t in ts)
write("constant/fvModels", f"""mixing
{{
    type            turbulentMixingDiffusion;
    phase           gas;
    species         (O2 CO);
    energy          {a.energy};
    tracer          N2;
    YAir            1;
    YProducts       0;
    zMin            1e-6;
    D
    {{
        type            table;
        values
        (
{table}
        );
    }}
}}
""")

# --------------------------------------------------------------------- 0
BC = """boundaryField
{{
    ends
    {{
        type            {t};
    }}
    sides
    {{
        type            empty;
    }}
}}
"""


def field(name, dims, vals, cls="volScalarField"):
    if isinstance(vals, str):
        internal = f"uniform {vals}"
    else:
        internal = ("nonuniform List<scalar>\n" + f"{len(vals)}\n(\n"
                    + "\n".join(f"{v:.12g}" for v in vals) + "\n)")
    bc = BC.format(t="symmetry")
    write(f"0/{name}", f"dimensions      {dims};\n\ninternalField   {internal};\n\n{bc}", cls)


left = [x < 0 for x in xc]
field("p", "[1 -1 -2 0 0 0 0]", f"{P0}")
T = [a.TLeft if l else a.TRight for l in left]
field("T", "[0 0 0 1 0 0 0]", T)
field("T.gas", "[0 0 0 1 0 0 0]", T)
field("U", "[0 1 -1 0 0 0 0]", "(0 0 0)", "volVectorField")
field("U.gas", "[0 1 -1 0 0 0 0]", "(0 0 0)", "volVectorField")
field("alpha.gas", "[0 0 0 0 0 0 0]", "1")
# Abel-Nobel: p = rho R T/(1 - b rho)  ->  rho = p/(R T + b p)
def rho(Y, T):
    """Abel-Nobel density of a N2/CO/O2 mixture (b, W of the thermo file)"""
    W = {"N2": 28.0134, "CO": 28.0101, "O2": 31.9988}
    b = {"N2": 1.3811e-3, "CO": 1.4116e-3, "O2": 9.9432e-4}
    R = sum(Y[s]*8314.47/W[s] for s in Y)
    return P0/(R*T + sum(Y[s]*b[s] for s in Y)*P0)


field("rho.gas", "[1 -3 0 0 0 0 0]",
      [rho({"N2": 1.0} if l else {"CO": 0.5, "O2": 0.5} if x > XO2 else {"CO": 1.0}, t)
       for l, x, t in zip(left, xc, T)])
field("Ydefault.gas", "[0 0 0 0 0 0 0]", "0")
# Mask check: for x > XO2 the products are half O2 (z = 1, outside the mixing
# layer, so the CO/O2 step there must stay sharp)
o2 = [0.5 if x > XO2 else 0.0 for x in xc]
field("O2.gas", "[0 0 0 0 0 0 0]", o2)
field("N2.gas", "[0 0 0 0 0 0 0]", [1.0 if l else 0.0 for l in left])
field("CO.gas", "[0 0 0 0 0 0 0]", [0.0 if l else 1.0 - y for l, y in zip(left, o2)])
