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

#include "eMeltingBlastThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::eMeltingThermo<EquationOfState>::eMeltingThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Cv_(dict.subDict("thermodynamics").lookup<scalar>("Cv")),
    CvLiquid_(dict.subDict("thermodynamics").lookup<scalar>("CvLiquid")),
    Tm_(dict.subDict("thermodynamics").lookup<scalar>("Tm")),
    Lm_(dict.subDict("thermodynamics").lookup<scalar>("Lm")),
    dTm_(dict.subDict("thermodynamics").lookupOrDefault<scalar>("dTm", 2.0)),
    Hf_(dict.subDict("thermodynamics").lookup<scalar>("Hf")),
    Tref_(dict.subDict("thermodynamics").lookupOrDefault<scalar>("Tref", Tstd)),
    Esref_
    (
        dict.subDict("thermodynamics").lookupOrDefault<scalar>
        (
            "Esref",
            Cv_*Tref_
        )
    )
{
    if (dTm_ <= 0 || Tref_ >= T1())
    {
        FatalIOErrorInFunction(dict)
            << "eMelting requires dTm > 0 and Tref < Tm - dTm/2"
            << " (dTm = " << dTm_ << ", Tref = " << Tref_
            << ", Tm = " << Tm_ << ")" << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::eMeltingThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Cv", Cv_);
    dict.add("CvLiquid", CvLiquid_);
    dict.add("Tm", Tm_);
    dict.add("Lm", Lm_);
    dict.add("dTm", dTm_);
    dict.add("Hf", Hf_);
    dict.add("Tref", Tref_);
    dict.add("Esref", Esref_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const eMeltingThermo<EquationOfState>& mt
)
{
    mt.write(os);
    return os;
}

// ************************************************************************* //
