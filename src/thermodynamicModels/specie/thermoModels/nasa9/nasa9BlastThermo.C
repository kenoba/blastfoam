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

\*---------------------------------------------------------------------------*/

#include "nasa9BlastThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
void Foam::nasa9Thermo<EquationOfState>::checkRanges
(
    const nasa9Thermo& nt1,
    const nasa9Thermo& nt2
)
{
    bool same = nt1.nRanges_ == nt2.nRanges_;
    if (same)
    {
        for (label i = 0; i <= nt1.nRanges_; i++)
        {
            same = same && !notEqual(nt1.Tswitch_[i], nt2.Tswitch_[i]);
        }
    }

    if (!same)
    {
        FatalErrorInFunction
            << "Temperature ranges " << nt1.Tswitch() << " for "
            << (nt1.name().size() ? nt1.name() : "others")
            << " != " << nt2.Tswitch() << " for "
            << (nt2.name().size() ? nt2.name() : "others") << nl
            << "    All species of a nasa9 mixture must use the same ranges"
            << exit(FatalError);
    }
}


template<class EquationOfState>
void Foam::nasa9Thermo<EquationOfState>::setRanges
(
    const UList<scalar>& Tswitch,
    const UList<coeffArray>& coeffs
)
{
    if (Tswitch.size() != coeffs.size() + 1)
    {
        FatalErrorInFunction
            << "Tswitch has " << Tswitch.size() << " entries for "
            << coeffs.size() << " coefficient ranges (expected "
            << coeffs.size() + 1 << ')'
            << exit(FatalError);
    }

    if (coeffs.size() < 1 || coeffs.size() > maxRanges_)
    {
        FatalErrorInFunction
            << coeffs.size() << " temperature ranges for "
            << this->name() << "; nasa9 supports 1 to " << maxRanges_
            << exit(FatalError);
    }

    for (label i = 0; i < Tswitch.size() - 1; i++)
    {
        if (Tswitch[i] >= Tswitch[i + 1])
        {
            FatalErrorInFunction
                << "Tswitch " << Tswitch << " is not increasing for "
                << this->name()
                << exit(FatalError);
        }
    }

    nRanges_ = coeffs.size();
    Tswitch_ = Zero;
    forAll(Tswitch, i)
    {
        Tswitch_[i] = Tswitch[i];
    }
    for (label i = 0; i < maxRanges_; i++)
    {
        coeffs_[i] = Zero;
    }
    forAll(coeffs, i)
    {
        coeffs_[i] = coeffs[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::nasa9Thermo<EquationOfState>::nasa9Thermo(const dictionary& dict)
:
    EquationOfState(dict),
    nRanges_(0)
{
    const dictionary& thermoDict = dict.subDict("thermodynamics");
    setRanges
    (
        List<scalar>(thermoDict.lookup("Tswitch")),
        List<coeffArray>(thermoDict.lookup("coeffs"))
    );

    // Convert coefficients to mass-basis
    for (label i = 0; i < nRanges_; i++)
    {
        for (label j = 0; j < nCoeffs_; j++)
        {
            coeffs_[i][j] *= this->R();
        }
    }
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const nasa9Thermo<EquationOfState>& nt
)
{
    nt.write(os);
    return os;
}


// ************************************************************************* //
