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

#include "AbelNobelSG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::AbelNobelSG<Specie>::AbelNobelSG
(
    const dictionary& dict
)
:
    Specie(dict),
    b_(dict.subDict("equationOfState").lookup<scalar>("b")),
    gamma_(dict.subDict("equationOfState").lookup<scalar>("gamma")),
    pInf_(dict.subDict("equationOfState").lookup<scalar>("pInf")),
    Cv_(dict.subDict("thermodynamics").lookup<scalar>("Cv")),
    pCav_
    (
        dict.subDict("equationOfState").lookupOrDefault<scalar>("pCav", 0.0)
    )
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::AbelNobelSG<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("b", b_);
    dict.add("gamma", gamma_);
    dict.add("pInf", pInf_);
    dict.add("pCav", pCav_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const AbelNobelSG<Specie>& an
)
{
    an.write(os);
    return os;
}


// ************************************************************************* //
