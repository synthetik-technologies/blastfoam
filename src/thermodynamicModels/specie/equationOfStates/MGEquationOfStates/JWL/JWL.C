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

#include "JWL.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::JWL<Specie>::JWL
(
    const dictionary& dict
)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    rhoCutOff_
    (
        dict.subDict("equationOfState").lookupOrDefault
        (
            "rhoCutOff",
            0.0
        )
    ),
    omega_(dict.subDict("equationOfState").lookup<scalar>("omega")),
    gammaIdeal_
    (
        dict.subDict("equationOfState").lookupOrDefault
        (
            "gammaIdeal",
            omega_ + 1.0
        )
    ),
    A_(dict.subDict("equationOfState").lookup<scalar>("A")),
    B_(dict.subDict("equationOfState").lookup<scalar>("B")),
    R1_(dict.subDict("equationOfState").lookup<scalar>("R1")),
    R2_(dict.subDict("equationOfState").lookup<scalar>("R2")),
    e0_(dict.subDict("equationOfState").lookupOrDefault<scalar>("e0", 0.0))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::JWL<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("rho0", rho0_);
    dict.add("omega", omega_);
    dict.add("A", A_);
    dict.add("B", B_);
    dict.add("R1", R1_);
    dict.add("R2", R2_);
    dict.add("e0", e0_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const JWL<Specie>& j
)
{
    j.write(os);
    return os;
}


// ************************************************************************* //
