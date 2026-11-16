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

Class
    Foam::fv::turbulentMixingDiffusion

Description
    Subgrid turbulent mixing of detonation products and air: a gradient
    diffusion term div(alpha rho D grad(phi)) in the species and energy
    equations of one fluid phase, restricted to the products/air mixing layer.

    The mixing layer is where the mixture fraction z of the products, from an
    inert tracer species (N2 by default),

        z = (Y_tracer - YAir)/(YProducts - YAir),

    lies between zMin and 1 - zMin. The face coefficient is interpolated
    from the masked cell values, so the layer can widen by diffusion. Outside
    the layer (shocks in air, the products core) nothing is added. D(t) is a
    Function1 of time [m^2/s] (e.g. from a mixing-layer growth law
    D = C w dw/dt, set by the case).

    The term is added in the implicit stage of blastEulerFoam (the species and
    energy equations of the phase call the fvModels with their field). It is
    in flux form, so species mass and energy are conserved; the energy is the
    phase (sensible) energy, i.e. unity Lewis number.

Usage
    \verbatim
    mixing
    {
        type            turbulentMixingDiffusion;
        phase           gas;
        species         (O2 N2 CO ...);
        energy          yes;
        tracer          N2;
        YAir            0.767;
        YProducts       0.183;
        zMin            1e-3;
        D               table ((0 0) (1e-3 0.5));
    }
    \endverbatim

SourceFiles
    turbulentMixingDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentMixingDiffusion_H
#define turbulentMixingDiffusion_H

#include "fvModel.H"
#include "Function1.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyTopoChangeMap;
class polyMeshMap;
class polyDistributionMap;

namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class turbulentMixingDiffusion Declaration
\*---------------------------------------------------------------------------*/

class turbulentMixingDiffusion
:
    public fvModel
{
    // Private Data

        //- Name of the phase
        word phaseName_;

        //- Diffused species (without the phase group)
        wordList species_;

        //- Diffuse the phase energy
        bool energy_;

        //- Name of the energy field (without the phase group)
        word energyName_;

        //- Tracer species defining the mixture fraction
        word tracer_;

        //- Tracer mass fraction in air and in the products
        scalar YAir_;
        scalar YProducts_;

        //- Mixture fraction bounds of the mixing layer
        scalar zMin_;

        //- Turbulent diffusivity [m^2/s] as a function of time
        autoPtr<Function1<scalar>> D_;

        //- Cached mixing-layer mask and the time index it belongs to
        mutable autoPtr<volScalarField> maskPtr_;
        mutable label maskTimeIndex_;


    // Private Member Functions

        //- Read the coefficients
        void readCoeffs();

        //- Mixing-layer mask (1 in the layer, 0 elsewhere), cached per step
        const volScalarField& mask() const;

        //- Face diffusion coefficient for the given (alpha) rho
        tmp<surfaceScalarField> Gamma(const volScalarField& alphaRho) const;

        //- Is the field one that is diffused
        bool diffuses(const word& fieldName) const;


public:

    //- Runtime type information
    TypeName("turbulentMixingDiffusion");


    // Constructors

        //- Construct from components
        turbulentMixingDiffusion
        (
            const word& name,
            const word& modelType,
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~turbulentMixingDiffusion() = default;


    // Member Functions

        //- Fields to which the model adds a source
        virtual wordList addSupFields() const;

        //- Does the model add a source to the field
        virtual bool addsSupToField(const word& fieldName) const;


        // Sources

            //- Add to a compressible equation (species: rho = alpha rho)
            virtual void addSup
            (
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const;

            //- Add to a phase equation (energy)
            virtual void addSup
            (
                const volScalarField& alpha,
                const volScalarField& rho,
                const volScalarField& field,
                fvMatrix<scalar>& eqn
            ) const;


        // Mesh motion and topology changes

            virtual bool movePoints()
            {
                maskPtr_.clear();
                return true;
            }

            virtual void topoChange(const polyTopoChangeMap&)
            {
                maskPtr_.clear();
            }

            virtual void mapMesh(const polyMeshMap&)
            {
                maskPtr_.clear();
            }

            virtual void distribute(const polyDistributionMap&)
            {
                maskPtr_.clear();
            }


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
