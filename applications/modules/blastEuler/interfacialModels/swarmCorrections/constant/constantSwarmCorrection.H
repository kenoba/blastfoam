/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2026-11-24 constant drag multiplier (blastfoam fork)
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::swarmCorrections::constantSwarm

Description
    Constant multiplier of the drag coefficient, e.g. for drag sensitivity
    studies or an empirical cloud (collective) drag factor:

    \verbatim
        swarmCorrection
        {
            type            constant;
            Cs              3;
        }
    \endverbatim

SourceFiles
    constantSwarmCorrection.C

\*---------------------------------------------------------------------------*/

#ifndef constantSwarmCorrection_H
#define constantSwarmCorrection_H

#include "swarmCorrection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace swarmCorrections
{

/*---------------------------------------------------------------------------*\
                           Class constantSwarm Declaration
\*---------------------------------------------------------------------------*/

class constantSwarm
:
    public swarmCorrection
{
    // Private data

        //- Drag multiplier
        const scalar Cs_;

public:

    //- Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from a dictionary and a phase pair
        constantSwarm
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~constantSwarm();


    // Member Functions

        //- Swarm correction coefficient
        virtual tmp<volScalarField> Cs() const;

        //- Swarm correction coefficient
        virtual scalar cellCs(const label celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace swarmCorrections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
