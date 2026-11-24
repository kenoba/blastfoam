/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2025-06-09 Jeff Heylmun:    Added cell based returns
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

\*---------------------------------------------------------------------------*/

#include "constantSwarmCorrection.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace swarmCorrections
{
    defineTypeNameAndDebug(constantSwarm, 0);
    addToRunTimeSelectionTable(swarmCorrection, constantSwarm, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swarmCorrections::constantSwarm::constantSwarm
(
    const dictionary& dict,
    const phasePair& pair
)
:
    swarmCorrection(dict, pair),
    Cs_(dict.lookup<scalar>("Cs"))
{
    if (Cs_ <= 0)
    {
        FatalIOErrorInFunction(dict)
            << "Cs must be positive, Cs = " << Cs_ << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swarmCorrections::constantSwarm::~constantSwarm()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::swarmCorrections::constantSwarm::Cs() const
{
    const fvMesh& mesh(this->pair_.phase1().mesh());

    return volScalarField::New
    (
        "Cs",
        mesh,
        dimensionedScalar("Cs", dimless, Cs_)
    );
}


Foam::scalar
Foam::swarmCorrections::constantSwarm::cellCs(const label) const
{
    return Cs_;
}


// ************************************************************************* //
