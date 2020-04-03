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

#include "solidJWL.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::solidJWL<Specie>::solidJWL
(
    const dictionary& dict
)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookupType<scalar>("rho0")),
    omega_(dict.subDict("equationOfState").lookupType<scalar>("omega")),
    A_(dict.subDict("equationOfState").lookupType<scalar>("A")),
    B_(dict.subDict("equationOfState").lookupType<scalar>("B")),
    R1_(dict.subDict("equationOfState").lookupType<scalar>("R1")),
    R2_(dict.subDict("equationOfState").lookupType<scalar>("R2"))
{}

// ************************************************************************* //
