/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
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

#include "reactingDiameterModel.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(reactingDiameterModel, 0);
    addToRunTimeSelectionTable(diameterModel, reactingDiameterModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::reactingDiameterModel::reactingDiameterModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    diameter(mesh, dict, phaseName),
    rate_(diameterReactionRate::New(*this, dict)),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    TName_(dict.lookupOrDefault("TName", IOobject::groupName("T", phaseName))),
    oxidiserPhase_(dict.lookupOrDefault("oxidiserPhase", word::null)),
    oxidiserSpecies_
    (
        dict.lookupOrDefault("oxidiserSpecies", wordList())
    ),
    YoxidiserRef_
    (
        oxidiserSpecies_.size()
      ? dict.lookup<scalar>("YoxidiserRef")
      : 1.0
    ),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 1e-4)),
    ignitedFraction_(dict.lookupOrDefault<bool>("ignitedFraction", false)),
    Tign_(ignitedFraction_ ? dict.lookup<scalar>("Tign") : great),
    tauIgn_
    (
        ignitedFraction_ ? dict.lookupOrDefault<scalar>("tauIgn", 0) : 0
    ),
    N_
    (
        IOobject
        (
            IOobject::groupName("N", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", inv(dimVolume), 1.0)
    ),
    d0_
    (
        IOobject
        (
            IOobject::groupName("d0", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->d_
    ),
    dVdt_
    (
        IOobject
        (
            IOobject::groupName("dVdt", phaseName),
            mesh.time().name(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimVolume/dimTime, 0.0)
    ),
    ignited_
    (
        IOobject
        (
            IOobject::groupName("ignited", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            ignitedFraction_ ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    ignitionProgress_
    (
        IOobject
        (
            IOobject::groupName("ignitionProgress", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            tauIgn_ > 0 ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless, 0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    if (tauIgn_ < 0)
    {
        FatalIOErrorInFunction(dict)
            << "tauIgn must not be negative" << exit(FatalIOError);
    }
    this->requireD();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::reactingDiameterModel::~reactingDiameterModel()
{}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::diameterModels::reactingDiameterModel::oxidiserFactor() const
{
    const fvMesh& mesh = this->d_.mesh();
    const volScalarField& alpha =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", this->d_.group())
        );
    const volScalarField& alphaOx =
        mesh.lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", oxidiserPhase_)
        );

    // Fraction of the fluid around the particles that carries the oxidiser
    tmp<volScalarField> tf
    (
        volScalarField::New
        (
            IOobject::groupName("oxidiserFactor", this->d_.group()),
            min
            (
                max(alphaOx, scalar(0))
               /max(1.0 - alpha, dimensionedScalar(dimless, small)),
                scalar(1)
            )
        )
    );

    if (oxidiserSpecies_.size())
    {
        volScalarField Yox
        (
            volScalarField::New
            (
                "Yoxidiser",
                mesh,
                dimensionedScalar(dimless, 0)
            )
        );
        forAll(oxidiserSpecies_, i)
        {
            Yox +=
                mesh.lookupObject<volScalarField>
                (
                    IOobject::groupName(oxidiserSpecies_[i], oxidiserPhase_)
                );
        }
        tf.ref() *= min(max(Yox, scalar(0))/YoxidiserRef_, scalar(1));
    }

    return tf;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::reactingDiameterModel::update()
{
    const volScalarField& alpha =
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", this->d_.group())
        );

    forAll(this->d_, celli)
    {
        if (alpha[celli] <= alphaMin_)
        {
            this->d_[celli] = d0_[celli];
        }
        else
        {
            this->d_[celli] = Foam::max(1e-10, Foam::min(this->d_[celli], d0_[celli]));
        }
    }
    this->d_.correctBoundaryConditions();

    const dimensionedScalar piBy6
    (
        "piBy6",
        dimless,
        constant::mathematical::pi/6.0
    );
    // Number density from the current diameter: the volume fraction
    // shrinks with the particles as they burn, so N = alpha/V(d0) would
    // make the mass loss dm/dt = m (dV/dt)/V0 (exponential, leaving 1/e of
    // the mass when d -> 0). The floor must be far below the volume of a
    // micron particle (5e-19 m^3), not small = 1e-15 m^3.
    N_ =
        alpha
       /max(piBy6*pow3(this->d_), dimensionedScalar(dimVolume, vSmall));
}


void Foam::diameterModels::reactingDiameterModel::solve()
{
    const volScalarField& p
    (
        this->d_.mesh().lookupObject<volScalarField>(pName_)
    );
    const volScalarField& T
    (
        this->d_.mesh().lookupObject<volScalarField>(TName_)
    );
    solve(p, T);
}


void Foam::diameterModels::reactingDiameterModel::solve
(
    const volScalarField& pi,
    const volScalarField& T
)
{
    // One stage of the multi-stage time integration. As for the conserved
    // fields of the phase system, the change is evaluated from the current
    // stage values, then
    //
    //     field = (blended old stage values) + dt (blended stage changes).
    //
    // Before 2026-11-21 the ignited fraction and the diameter were advanced
    // from the unblended current stage value, which with RK2SSP transported
    // them 1.5 times faster than the particles (the ignited fraction left
    // the particle cloud and was reset), and the mass transfer included the
    // fall of the cell-mean diameter by transport, with the diameter not
    // allowed to rise, which removed particles numerically (1D FOI shot 5:
    // 23% of the particles by 1 ms, 59% of the Al mass by 16 ms, of which
    // shrinking explained 9%).

    const volScalarField& alpha =
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("alpha", this->d_.group())
        );
    const surfaceScalarField& phi
    (
        this->d_.mesh().lookupObject<surfaceScalarField>
        (
            IOobject::groupName("phi", this->d_.group())
        )
    );
    const dimensionedScalar& dT(this->d_.time().deltaT());
    const word divScheme("div(" + phi.name() + ',' + this->d_.name() + ')');

    // Burn rate of the particles at the current stage
    volScalarField dDdt(-2.0*rate_->dDdt(pi, T));
    if (!oxidiserPhase_.empty())
    {
        dDdt *= oxidiserFactor();
    }
    if (ignitedFraction_)
    {
        dDdt *= ignited_;
    }
    dDdt.min(0.0);

    // Mass transfer: volume change of a particle by the burning alone,
    // dV/dt = (pi/2) d^2 dd/dt, at most the particle volume in one step
    const dimensionedScalar piBy6(dimless, constant::mathematical::pi/6.0);
    dVdt_ = 3.0*piBy6*sqr(this->d_)*dDdt;
    dVdt_ = max(dVdt_, -piBy6*pow3(this->d_)/dT);

    // Ignited fraction: transported with the particles, then set where the
    // particles have reached Tign (tauIgn = 0) or have been at or above Tign
    // for tauIgn (ignition progress, transported with the particles)
    if (ignitedFraction_)
    {
        const bool delay = tauIgn_ > 0;
        if (delay)
        {
            volScalarField deltaProgress
            (
                IOobject::groupName("ignitionProgressDelta", this->d_.group()),
                fvc::div(phi, ignitionProgress_, divScheme)
              - ignitionProgress_*fvc::div(phi)
              - pos0(T - dimensionedScalar(dimTemperature, Tign_))
               /dimensionedScalar(dimTime, tauIgn_)
            );
            this->storeAndBlendOld(ignitionProgress_);
            this->storeAndBlendDelta(deltaProgress);
            ignitionProgress_ -= dT*deltaProgress;
            forAll(ignitionProgress_, celli)
            {
                ignitionProgress_[celli] =
                    alpha[celli] <= alphaMin_
                  ? 0
                  : Foam::min(Foam::max(ignitionProgress_[celli], 0.0), 1.0);
            }
            ignitionProgress_.correctBoundaryConditions();
        }

        volScalarField deltaPhiIgn
        (
            IOobject::groupName("ignitedDelta", this->d_.group()),
            fvc::div(phi, ignited_, divScheme) - ignited_*fvc::div(phi)
        );
        this->storeAndBlendOld(ignited_);
        this->storeAndBlendDelta(deltaPhiIgn);
        ignited_ -= dT*deltaPhiIgn;
        forAll(ignited_, celli)
        {
            if (alpha[celli] <= alphaMin_)
            {
                ignited_[celli] = 0;
            }
            else if
            (
                delay
              ? ignitionProgress_[celli] >= 1 - small
              : T[celli] >= Tign_
            )
            {
                ignited_[celli] = 1;
            }
            else
            {
                ignited_[celli] = Foam::min(Foam::max(ignited_[celli], 0.0), 1.0);
            }
        }
        ignited_.correctBoundaryConditions();
    }

    // Diameter: transported with the particles (it may rise where larger
    // particles are carried in) and reduced by the burning
    volScalarField ddt
    (
        IOobject::groupName("dDelta", this->d_.group()),
        dDdt - (fvc::div(phi, this->d_, divScheme) - this->d_*fvc::div(phi))
    );
    this->storeAndBlendOld(this->d_);
    this->storeAndBlendDelta(ddt);
    this->d_ += dT*ddt;
    forAll(this->d_, celli)
    {
        if (alpha[celli] <= alphaMin_)
        {
            this->d_[celli] = d0_[celli];
        }
        else
        {
            this->d_[celli] = Foam::max(1e-10, Foam::min(this->d_[celli], d0_[celli]));
        }
    }
    this->d_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::reactingDiameterModel::dMdt() const
{
    const volScalarField& rho
    (
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("rho", this->d_.group())
        )
    );

    return dVdt_*rho*N_;
}

// ************************************************************************* //
