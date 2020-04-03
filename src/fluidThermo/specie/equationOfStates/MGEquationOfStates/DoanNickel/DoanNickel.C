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

#include "DoanNickel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::DoanNickel<Specie>::DoanNickel
(
    const dictionary& dict
)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookupType<scalar>("rho0")),
    e0_(dict.subDict("equationOfState").lookupType<scalar>("e0")),
    a_(dict.subDict("equationOfState").lookupType<scalar>("a")),
    b_(dict.subDict("equationOfState").lookupType<scalar>("b")),
    c_(dict.subDict("equationOfState").lookupType<scalar>("c")),
    d_(dict.subDict("equationOfState").lookupType<scalar>("d")),
    g_(dict.subDict("equationOfState").lookupType<scalar>("g")),

    a1_(dict.subDict("equationOfState").lookupType<scalar>("a1")),
    a2_(dict.subDict("equationOfState").lookupType<scalar>("a2")),
    a3_(dict.subDict("equationOfState").lookupType<scalar>("a3")),

    E1_(dict.subDict("equationOfState").lookupType<scalar>("E1")),
    E11_(dict.subDict("equationOfState").lookupType<scalar>("E11")),
    E1Offset_(dict.subDict("equationOfState").lookupType<scalar>("E1Offset")),
    deltaE1_(dict.subDict("equationOfState").lookupType<scalar>("deltaE1")),
    deltaE1Pow_(dict.subDict("equationOfState").lookupType<scalar>("deltaE1Pow")),

    E2_(dict.subDict("equationOfState").lookupType<scalar>("E2")),
    E22_(dict.subDict("equationOfState").lookupType<scalar>("E22")),
    deltaE2_(dict.subDict("equationOfState").lookupType<scalar>("deltaE2")),
    E2Pow_(dict.subDict("equationOfState").lookupType<scalar>("E2Pow")),
    deltaE2Pow_(dict.subDict("equationOfState").lookupType<scalar>("deltaE2Pow")),

    E3_(dict.subDict("equationOfState").lookupType<scalar>("E3")),
    E33_(dict.subDict("equationOfState").lookupType<scalar>("E33")),
    deltaE3_(dict.subDict("equationOfState").lookupType<scalar>("deltaE3"))
{}

// ************************************************************************* //
