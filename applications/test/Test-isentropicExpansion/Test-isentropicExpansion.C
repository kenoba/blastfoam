/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2026
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-isentropicExpansion

Description
    0D closed-volume isentropic expansion of a nasa9/Abel-Nobel gas mixture
    (fixed composition) with the compiled thermo classes, as the flow solver
    evaluates them: the first law de = -p dv is integrated (trapezoidal rule,
    log-spaced volumes) from (rho0, T0) to rho1, and the temperature at
    each step follows from the energy with TRhoE. Prints the mixture gas
    constant and covolume, the species Cp, Cv, Es and pressure at a test
    state, and the expansion (rho T p e Cv).

    Test-isentropicExpansion -file thermo.nasa9Gas \
        -mixture '((N2 0.18) (CO 0.27) ("C(gr)" 0.18) ...)' \
        -rho0 2000 -T0 3600 -rho1 0.1 [-nSteps 20000]

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Tuple2.H"
#include "specieBlast.H"
#include "AbelNobel.H"
#include "nasa9BlastThermo.H"
#include "thermoModel.H"

using namespace Foam;

typedef thermoModel<nasa9Thermo<AbelNobel<specieBlast>>> ThermoType;


// Species dictionary with a tight temperature tolerance
dictionary specieDict(const dictionary& dict, const word& name)
{
    dictionary d(dict.subDict(name));
    d.subDict("thermodynamics").set("relTol", 1e-12);
    return d;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("file", "species file");
    argList::addOption("mixture", "list of (species massFraction)");
    argList::addOption("rho0", "initial density [kg/m^3]");
    argList::addOption("T0", "initial temperature [K]");
    argList::addOption("rho1", "final density [kg/m^3]");
    argList::addOption("nSteps", "number of volume steps (default 20000)");
    argList args(argc, argv, false, true);
    IOstream::defaultPrecision(14);
    Sout.precision(14);
    Info().precision(14);

    IFstream is(args.optionRead<fileName>("file"));
    const dictionary dict(is);
    const List<Tuple2<word, scalar>> Y
    (
        args.optionRead<List<Tuple2<word, scalar>>>("mixture")
    );
    const scalar rho0 = args.optionRead<scalar>("rho0");
    const scalar T0 = args.optionRead<scalar>("T0");
    const scalar rho1 = args.optionRead<scalar>("rho1");
    const label nSteps = args.optionLookupOrDefault<label>("nSteps", 20000);

    // Species at a test state
    const scalar rhoTest = 10, TTest = 1500;
    Info<< "# specie name Rg b Cp Cv Es p (rho " << rhoTest << ", T "
        << TTest << ")" << endl;
    autoPtr<ThermoType> mix;
    forAll(Y, i)
    {
        const ThermoType t(specieDict(dict, Y[i].first()));
        Info<< "specie " << Y[i].first() << ' ' << t.Rg() << ' ' << t.b()
            << ' ' << t.Cp(rhoTest, 0, TTest) << ' '
            << t.Cv(rhoTest, 0, TTest) << ' ' << t.Es(rhoTest, 0, TTest)
            << ' ' << t.p(rhoTest, 0, TTest, false) << endl;
        if (mix.valid())
        {
            mix() += Y[i].second()*t;
        }
        else
        {
            mix.set(new ThermoType("mixture", Y[i].second()*t));
        }
    }
    const ThermoType& t = mix();
    Info<< "mixture " << t.Rg() << ' ' << t.b() << ' '
        << t.Cp(rhoTest, 0, TTest) << ' ' << t.Cv(rhoTest, 0, TTest) << ' '
        << t.Es(rhoTest, 0, TTest) << ' ' << t.p(rhoTest, 0, TTest, false)
        << endl;

    // Expansion: e_{n+1} = e_n - (p_n + p_{n+1})/2 (v_{n+1} - v_n)
    Info<< "# state rho T p e Cv" << endl;
    scalar rho = rho0;
    scalar T = T0;
    scalar e = t.Es(rho, 0, T);
    scalar p = t.p(rho, e, T, false);
    Info<< "state " << rho << ' ' << T << ' ' << p << ' ' << e << ' '
        << t.Cv(rho, e, T) << endl;
    const scalar r = std::pow(rho1/rho0, scalar(1)/nSteps);
    const label nPrint = max(nSteps/200, 1);
    for (label n = 1; n <= nSteps; n++)
    {
        const scalar rhoNew = rho0*Foam::pow(r, scalar(n));
        const scalar dv = 1/rhoNew - 1/rho;
        scalar TNew = T;
        scalar eNew = e - p*dv;
        for (label iter = 0; iter < 50; iter++)
        {
            TNew = t.TRhoE(TNew, rhoNew, eNew);
            const scalar pNew = t.p(rhoNew, eNew, TNew, false);
            const scalar eNext = e - 0.5*(p + pNew)*dv;
            if (mag(eNext - eNew) < 1e-12*max(mag(e), 1))
            {
                eNew = eNext;
                break;
            }
            eNew = eNext;
        }
        TNew = t.TRhoE(TNew, rhoNew, eNew);
        rho = rhoNew;
        T = TNew;
        e = eNew;
        p = t.p(rho, e, T, false);
        if (n % nPrint == 0 || n == nSteps)
        {
            Info<< "state " << rho << ' ' << T << ' ' << p << ' ' << e << ' '
                << t.Cv(rho, e, T) << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
