/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "eConstBlastThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::eConstThermo<EquationOfState>::eConstThermo(const dictionary& dict)
:
    EquationOfState(dict),
    Cv_(dict.subDict("thermodynamics").lookup<scalar>("Cv")),
    Hf_(dict.subDict("thermodynamics").lookup<scalar>("Hf")),
    Tref_(dict.subDict("thermodynamics").lookupOrDefault<scalar>("Tref", 0)),
    Esref_(dict.subDict("thermodynamics").lookupOrDefault<scalar>("Esref", 0)),
    flameT_(dict.subDict("thermodynamics").lookupOrDefault("flameT", Hf_/Cv_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
void Foam::eConstThermo<EquationOfState>::write(Ostream& os) const
{
    EquationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("Cv", Cv_);
    dict.add("Hf", Hf_);
    if (Tref_ != Tstd)
    {
        dict.add("Tref", Tref_);
    }
    if (Esref_ != 0)
    {
        dict.add("Esref", Esref_);
    }
    if (flameT_ != 0)
    {
        dict.add("flameT", flameT_);
    }
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const eConstThermo<EquationOfState>& et
)
{
    et.write(os);
    return os;
}

// ************************************************************************* //
