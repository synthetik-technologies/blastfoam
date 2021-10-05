/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "BirchMurnaghan2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::BirchMurnaghan2<Specie>::BirchMurnaghan2
(
    const dictionary& dict
)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    pRef_(dict.subDict("equationOfState").lookup<scalar>("pRef")),
    K0_(dict.subDict("equationOfState").lookup<scalar>("K0")),
    Gamma_(dict.subDict("equationOfState").lookup<scalar>("Gamma"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::BirchMurnaghan2<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("rho0", rho0_);
    dict.add("pRef", pRef_);
    dict.add("K0", K0_);
    dict.add("K0Prime", K0Prime_);
    dict.add("Gamma", Gamma_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const BirchMurnaghan2<Specie>& bm2
)
{
    bm2.write(os);
    return os;
}


// ************************************************************************* //
