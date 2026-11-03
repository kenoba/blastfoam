/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "PaSRBlastCombustionModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
#include "wedgeFvPatch.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(PaSR, 0);
    addToRunTimeSelectionTable(blastCombustionModel, PaSR, fluid);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::PaSR::PaSR
(
    const word& modelType,
    const multicomponentBlastThermo& thermo,
    const bool isFluid,
    const word& combustionProperties
)
:
    laminar(modelType, thermo, isFluid, combustionProperties),
    Cmix_(this->coeffs().template lookup<scalar>("Cmix")),
    deltaMixing_(false),
    Cs_(0.17),
    kappa_
    (
        IOobject
        (
            thermo.phasePropertyName(typedName("kappa")),
            this->mesh().time().name(),
            this->mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh(),
        dimensionedScalar(dimless, 0)
    )
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::PaSR::~PaSR()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::combustionModels::PaSR::readCoeffs()
{
    const word mixingTime
    (
        this->coeffs().lookupOrDefault<word>("mixingTime", "turbulent")
    );

    if (mixingTime == "delta")
    {
        deltaMixing_ = true;
    }
    else if (mixingTime == "turbulent")
    {
        deltaMixing_ = false;
    }
    else
    {
        FatalIOErrorInFunction(this->coeffs())
            << "Unknown mixingTime " << mixingTime
            << ", valid options are turbulent and delta"
            << exit(FatalIOError);
    }

    Cs_ = this->coeffs().lookupOrDefault<scalar>("Cs", 0.17);
}


Foam::tmp<Foam::scalarField> Foam::combustionModels::PaSR::delta() const
{
    const fvMesh& mesh = this->mesh();
    const vectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceScalarField& magSf = mesh.magSf();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

    scalarField sumD(mesh.nCells(), 0);
    scalarField nD(mesh.nCells(), 0);

    forAll(own, facei)
    {
        const vector n(Sf[facei]/magSf[facei]);
        sumD[own[facei]] += mag((Cf[facei] - C[own[facei]]) & n);
        sumD[nei[facei]] += mag((Cf[facei] - C[nei[facei]]) & n);
        nD[own[facei]] += 1;
        nD[nei[facei]] += 1;
    }

    forAll(mesh.boundary(), patchi)
    {
        const fvPatch& p = mesh.boundary()[patchi];
        if (isA<wedgeFvPatch>(p) || isA<emptyFvPatch>(p))
        {
            continue;
        }

        const labelUList& fc = p.faceCells();
        const vectorField& pCf = p.Cf();
        const vectorField pn(p.nf());
        forAll(fc, i)
        {
            sumD[fc[i]] += mag((pCf[i] - C[fc[i]]) & pn[i]);
            nD[fc[i]] += 1;
        }
    }

    return 2*sumD/max(nD, scalar(1));
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::PaSR::correct()
{
    laminar::correct();

    if (deltaMixing_)
    {
        const compressibleMomentumTransportModel& turb = this->turbulence();

        tmp<volScalarField> tnuEff(turb.nuEff());
        const scalarField& nuEff = tnuEff();

        const volTensorField gradU(fvc::grad(turb.U()));
        const scalarField magD
        (
            sqrt(2.0)*mag(dev(symm(gradU.primitiveField())))
        );

        tmp<volScalarField> ttc(this->chemistryPtr_->tc());
        const scalarField& tc = ttc();

        const scalarField d(delta());

        forAll(kappa_, i)
        {
            const scalar nuSgs = sqr(Cs_*d[i])*magD[i];
            const scalar tk = Cmix_*sqr(d[i])/(nuEff[i] + nuSgs + vSmall);
            kappa_[i] = tk > small ? tc[i]/(tc[i] + tk) : 1.0;
        }

        return;
    }

    tmp<volScalarField> tepsilon(this->turbulence().epsilon());
    const scalarField& epsilon = tepsilon();

    tmp<volScalarField> tnuEff(this->turbulence().nuEff());
    const scalarField& nuEff = tnuEff();

    tmp<volScalarField> ttc(this->chemistryPtr_->tc());
    const scalarField& tc = ttc();

    forAll(epsilon, i)
    {
        const scalar tk =
            Cmix_*sqrt(max(nuEff[i]/(epsilon[i] + small), 0));

        if (tk > small)
        {
            kappa_[i] = tc[i]/(tc[i] + tk);
        }
        else
        {
            kappa_[i] = 1.0;
        }
    }
}


Foam::tmp<Foam::volScalarField::Internal>
Foam::combustionModels::PaSR::R(const label speciei) const
{
    return kappa_*laminar::R(speciei);
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::PaSR::R(volScalarField& Y) const
{
    return kappa_*laminar::R(Y);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::PaSR::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typedName("Qdot")),
        kappa_*laminar::Qdot()
    );
}


bool Foam::combustionModels::PaSR::read()
{
    if (laminar::read())
    {
        this->coeffs().lookup("Cmix") >> Cmix_;
        readCoeffs();
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
