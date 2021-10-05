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

#include "CochranChan.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::CochranChan<Specie>::CochranChan(const dictionary& dict)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    Gamma0_(dict.subDict("equationOfState").lookup<scalar>("Gamma0")),
    A_(dict.subDict("equationOfState").lookup<scalar>("A")),
    Epsilon1_(dict.subDict("equationOfState").lookup<scalar>("Epsilon1")),
    B_(dict.subDict("equationOfState").lookup<scalar>("B")),
    Epsilon2_(dict.subDict("equationOfState").lookup<scalar>("Epsilon2")),
    e0_(dict.subDict("equationOfState").lookupOrDefault<scalar>("e0", 0.0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::CochranChan<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("rho0", rho0_);
    dict.add("Gamma0", Gamma0_);
    dict.add("A", A_);
    dict.add("Epsilon1", Epsilon1_);
    dict.add("B", B_);
    dict.add("Epsilon2", Epsilon2_);
    dict.add("e0", e0_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CochranChan<Specie>& cc
)
{
    cc.write(os);
    return os;
}


// ************************************************************************* //
