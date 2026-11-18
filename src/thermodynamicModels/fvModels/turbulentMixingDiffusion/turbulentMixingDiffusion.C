/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
  \\    /   O peration     | Website:  https://openfoam.org
   \\  /    A nd           | Copyright (C) 2026
    \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "turbulentMixingDiffusion.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "fvmLaplacian.H"
#include "surfaceInterpolate.H"
#include "Time.H"
#include "fluidBlastThermo.H"
#include "basicSpecieBlastMixture.H"
#include "physicalProperties.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(turbulentMixingDiffusion, 0);
    addToRunTimeSelectionTable(fvModel, turbulentMixingDiffusion, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::turbulentMixingDiffusion::readCoeffs()
{
    phaseName_ = coeffs().lookup<word>("phase");
    species_ = coeffs().lookup<wordList>("species");
    energy_ = coeffs().lookupOrDefault<Switch>("energy", true);
    energyName_ = coeffs().lookupOrDefault<word>("energyName", "e");
    tracer_ = coeffs().lookupOrDefault<word>("tracer", "N2");
    YAir_ = coeffs().lookup<scalar>("YAir");
    YProducts_ = coeffs().lookup<scalar>("YProducts");
    zMin_ = coeffs().lookupOrDefault<scalar>("zMin", 1e-3);

    if (mag(YProducts_ - YAir_) < small)
    {
        FatalIOErrorInFunction(coeffs())
            << "YAir and YProducts of the tracer " << tracer_
            << " must differ" << exit(FatalIOError);
    }

    D_.reset
    (
        Function1<scalar>::New
        (
            "D",
            mesh().time().userUnits(),
            dimArea/dimTime,
            coeffs()
        ).ptr()
    );

    maskPtr_.clear();
    maskTimeIndex_ = -1;
    speciesChecked_ = false;
}


void Foam::fv::turbulentMixingDiffusion::checkSpecies() const
{
    if (speciesChecked_)
    {
        return;
    }
    speciesChecked_ = true;

    const word thermoName
    (
        IOobject::groupName(physicalProperties::typeName, phaseName_)
    );
    if (!mesh().foundObject<fluidBlastThermo>(thermoName))
    {
        WarningInFunction
            << "No thermo " << thermoName << " found; the species list of "
            << name() << " is not checked" << endl;
        return;
    }

    const basicSpecieBlastMixture* mixturePtr =
        dynamic_cast<const basicSpecieBlastMixture*>
        (
            &mesh().lookupObject<fluidBlastThermo>(thermoName)
        );
    if (!mixturePtr)
    {
        return;
    }

    wordList missing;
    forAll(mixturePtr->species(), i)
    {
        const word& s = mixturePtr->species()[i];
        if (findIndex(species_, s) == -1)
        {
            missing.append(s);
        }
    }

    if (missing.size())
    {
        FatalIOErrorInFunction(coeffs())
            << "Species " << missing << " of phase " << phaseName_
            << " are not in the species list of " << name() << nl
            << "    All species, the default specie included, must be "
            << "diffused: the thermo normalises the mass fractions when the "
            << "phase receives mass transfer, which converts an undiffused "
            << "species into the others." << exit(FatalIOError);
    }
}


const Foam::volScalarField& Foam::fv::turbulentMixingDiffusion::mask() const
{
    const label timeIndex = mesh().time().timeIndex();

    if
    (
        !maskPtr_.valid()
     || maskTimeIndex_ != timeIndex
     || maskPtr_->size() != mesh().nCells()
    )
    {
        const volScalarField& Yt =
            mesh().lookupObject<volScalarField>
            (
                IOobject::groupName(tracer_, phaseName_)
            );

        const volScalarField z((Yt - YAir_)/(YProducts_ - YAir_));

        maskPtr_.reset
        (
            volScalarField::New
            (
                IOobject::groupName("mixingMask", phaseName_),
                pos(z - zMin_)*pos(1.0 - zMin_ - z)
            ).ptr()
        );
        maskTimeIndex_ = timeIndex;
    }

    return maskPtr_();
}


Foam::tmp<Foam::surfaceScalarField>
Foam::fv::turbulentMixingDiffusion::Gamma
(
    const volScalarField& alphaRho
) const
{
    const dimensionedScalar D
    (
        dimArea/dimTime,
        max(D_->value(mesh().time().value()), 0)
    );

    return fvc::interpolate(alphaRho*mask())*D;
}


bool Foam::fv::turbulentMixingDiffusion::diffuses
(
    const word& fieldName
) const
{
    if (IOobject::group(fieldName) != phaseName_)
    {
        return false;
    }

    const word member(IOobject::member(fieldName));

    if (energy_ && member == energyName_)
    {
        return true;
    }

    return findIndex(species_, member) != -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::turbulentMixingDiffusion::turbulentMixingDiffusion
(
    const word& name,
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    fvModel(name, modelType, mesh, dict),
    phaseName_(),
    species_(),
    energy_(true),
    energyName_("e"),
    tracer_("N2"),
    YAir_(0),
    YProducts_(1),
    zMin_(1e-3),
    D_(),
    maskPtr_(),
    maskTimeIndex_(-1),
    speciesChecked_(false)
{
    readCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::wordList Foam::fv::turbulentMixingDiffusion::addSupFields() const
{
    wordList fields(species_.size() + (energy_ ? 1 : 0));

    forAll(species_, i)
    {
        fields[i] = IOobject::groupName(species_[i], phaseName_);
    }

    if (energy_)
    {
        fields.last() = IOobject::groupName(energyName_, phaseName_);
    }

    return fields;
}


bool Foam::fv::turbulentMixingDiffusion::addsSupToField
(
    const word& fieldName
) const
{
    return diffuses(fieldName);
}


void Foam::fv::turbulentMixingDiffusion::addSup
(
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (!diffuses(field.name()))
    {
        return;
    }

    checkSpecies();

    eqn += fvm::laplacian(Gamma(rho), field);
}


void Foam::fv::turbulentMixingDiffusion::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const volScalarField& field,
    fvMatrix<scalar>& eqn
) const
{
    if (!diffuses(field.name()))
    {
        return;
    }

    checkSpecies();

    eqn += fvm::laplacian(Gamma(alpha*rho), field);
}


bool Foam::fv::turbulentMixingDiffusion::read(const dictionary& dict)
{
    if (fvModel::read(dict))
    {
        readCoeffs();
        return true;
    }

    return false;
}


// ************************************************************************* //
