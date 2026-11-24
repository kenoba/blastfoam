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

Class
    Foam::diameterModels::reactingDiameterModel

Description
    Reacting diameter model used for mass transfer

    With an optional oxidiser phase the burn rate is scaled by the fraction
    of the surrounding fluid volume that this phase occupies and, with
    oxidiser species, by min(1, sum of their mass fractions/YoxidiserRef):

    \verbatim
        oxidiserPhase   gas;
        oxidiserSpecies (O2 CO2 H2O);   // optional
        YoxidiserRef    0.24;           // full rate at this mass fraction
    \endverbatim

    Without oxidiserPhase the particles burn wherever the rate model allows.

    With ignitedFraction the ignition is tracked by the field ignited.<phase>,
    the fraction of the particles that have ignited. It is transported with
    the particle volume flux (as the diameter), set to one where the particle
    temperature reaches Tign, and reset to zero below alphaMin; the burn rate
    of the rate model (applied to burning particles, without its own
    ignition check) is multiplied by it. Without it, a rate model such as
    Beckstead keeps burning wherever the cell-mean diameter is below d0, so
    the burning spreads with the mixing of the diameter field.

    With tauIgn > 0 the particles ignite only after being at or above Tign
    for tauIgn: the ignition progress ignitionProgress.<phase> grows at
    1/tauIgn where T >= Tign, is transported with the particles like the
    ignited fraction (so it accumulates along the particle paths and is kept
    when they cool), and the ignited fraction is set to one where it reaches
    one. This delays the ignition of particles heated in the early detonation
    products to the time they are in the fireball.

    The delay may depend on a temperature T' and on the particle diameter,

        tau = tauIgn exp(TaIgn/T') (d/dRefIgn)^nIgn,

    with T' the particle temperature or the field TIgnDelay. With
    TIgnDelay T.gas, Tign 0, tauIgn 2.5e-11 s and TaIgn 30000 K this is the
    ignition model of Kuhl and Boiko (LLNL-CONF-427973, 2010; Balakrishnan
    et al., Shock Waves 22, 2012), a fit to shock-tube delays of 3-5 um Al
    spheres at 1670-2170 K (0.43 ms at 1800 K); nIgn 2 with dRefIgn 4e-6
    scales it to other sizes (the heating time scales as d^2, a thermal
    explosion model gives d^2.3). That delay includes the heating of the
    particles, which the model resolves, so with the particle temperature
    T' and Tign near the melting point it counts twice.

    \verbatim
        ignitedFraction true;       // optional, default false
        Tign            2050;       // ignition temperature [K]
        tauIgn          5e-4;       // optional ignition delay [s], default 0
        TaIgn           0;          // optional activation temperature [K]
        nIgn            0;          // optional diameter exponent
        dRefIgn         4e-6;       // reference diameter [m] (nIgn != 0)
        TIgnDelay       T.gas;      // optional, default the particle T
    \endverbatim

SourceFiles
    reactingDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef reactingDiameterModel_H
#define reactingDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "diameterModel.H"
#include "diameterReactionRate.H"

namespace Foam
{
namespace diameterModels
{

/*---------------------------------------------------------------------------*\
                           Class reactingDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class reactingDiameterModel
:
    public diameter
{
// Private data

    //- Reaction rate model
    autoPtr<diameterReactionRate> rate_;

    //- Name of pressure field
    const word pName_;

    //- Name of temperature field
    const word TName_;

    //- Phase carrying the oxidiser (empty: no oxidiser limit)
    const word oxidiserPhase_;

    //- Oxidiser species in the oxidiser phase (empty: the whole phase)
    const wordList oxidiserSpecies_;

    //- Oxidiser mass fraction that gives the full burn rate
    const scalar YoxidiserRef_;

    //- Volume fraction below which the diameter is reset to d0 and the
    //  particles do not burn
    const scalar alphaMin_;

    //- Track the ignited fraction of the particles
    const bool ignitedFraction_;

    //- Ignition temperature (with ignitedFraction)
    const scalar Tign_;

    //- Ignition delay: time at or above Tign before ignition [s]
    //  (pre-exponential factor with TaIgn)
    const scalar tauIgn_;

    //- Activation temperature of the ignition delay [K]
    const scalar TaIgn_;

    //- Diameter exponent of the ignition delay
    const scalar nIgn_;

    //- Reference diameter of the ignition delay [m]
    const scalar dRefIgn_;

    //- Temperature field of the ignition delay
    const word TIgnDelayName_;

    //- Number density
    volScalarField N_;

    //- Initial diameter field
    volScalarField d0_;

    //- Stored change in volume per unit time
    volScalarField dVdt_;

    //- Fraction of the particles that have ignited
    volScalarField ignited_;

    //- Ignition progress (time at or above Tign over tauIgn)
    volScalarField ignitionProgress_;


    // Private Member Functions

        //- Burn rate scaling by the available oxidiser, in [0, 1]
        tmp<volScalarField> oxidiserFactor() const;

public:

    //- Runtime type information
    TypeName("reacting");

    // Constructor
    reactingDiameterModel
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );


    //- Destructor
    virtual ~reactingDiameterModel();

    // ODE functions

        //- Update number density
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Solve sub-step stepi given pressure and temperature fields
        virtual void solve(const volScalarField& p, const volScalarField& T);

        //- Return initial unburned diameter field
        virtual const volScalarField& d0() const
        {
            return d0_;
        }

        //- Return initial unburned diameter for a cell
        virtual scalar celld0(const label celli) const
        {
            return d0_[celli];
        }

        //- Return change in radius per unit time
        //  Only used in mass transfer
        virtual tmp<volScalarField> dMdt() const;

        //- Reaction rate
        virtual const diameterReactionRate* reactionRate() const
        {
            return rate_.valid() ? &rate_() : nullptr;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
