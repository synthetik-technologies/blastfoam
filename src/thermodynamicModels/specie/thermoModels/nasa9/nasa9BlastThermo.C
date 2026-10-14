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
    bool same = nt1.Tswitch_.size() == nt2.Tswitch_.size();
    if (same)
    {
        forAll(nt1.Tswitch_, i)
        {
            same = same && !notEqual(nt1.Tswitch_[i], nt2.Tswitch_[i]);
        }
    }

    if (!same)
    {
        FatalErrorInFunction
            << "Temperature ranges " << nt1.Tswitch_ << " for "
            << (nt1.name().size() ? nt1.name() : "others")
            << " != " << nt2.Tswitch_ << " for "
            << (nt2.name().size() ? nt2.name() : "others") << nl
            << "    All species of a nasa9 mixture must use the same ranges"
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::nasa9Thermo<EquationOfState>::nasa9Thermo(const dictionary& dict)
:
    EquationOfState(dict),
    Tswitch_(dict.subDict("thermodynamics").lookup("Tswitch")),
    coeffs_(dict.subDict("thermodynamics").lookup("coeffs"))
{
    if (Tswitch_.size() != coeffs_.size() + 1)
    {
        FatalIOErrorInFunction(dict)
            << "Tswitch has " << Tswitch_.size() << " entries for "
            << coeffs_.size() << " coefficient ranges (expected "
            << coeffs_.size() + 1 << ')'
            << exit(FatalIOError);
    }

    for (label i = 0; i < Tswitch_.size() - 1; i++)
    {
        if (Tswitch_[i] >= Tswitch_[i + 1])
        {
            FatalIOErrorInFunction(dict)
                << "Tswitch " << Tswitch_ << " is not increasing"
                << exit(FatalIOError);
        }
    }

    // Convert coefficients to mass-basis
    forAll(coeffs_, i)
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
