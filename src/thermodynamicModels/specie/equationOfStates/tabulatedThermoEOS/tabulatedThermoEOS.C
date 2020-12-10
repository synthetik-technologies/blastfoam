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

#include "tabulatedThermoEOS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::tabulatedThermoEOS<Specie>::tabulatedThermoEOS
(
    const dictionary& dict
)
:
    Specie(dict),
    pTable_
    (
        dict.subDict("equationOfState").lookup<fileName>("file"),
        dict.subDict("equationOfState").lookup<word>("mod"),
        dict.subDict("equationOfState").lookup<word>("rhoMod"),
        dict.subDict("equationOfState").lookup<word>("eMod"),
        dict.subDict("equationOfState").lookup<label>("nRho"),
        dict.subDict("equationOfState").lookup<label>("ne"),
        dict.subDict("equationOfState").lookup<scalar>("minRho"),
        dict.subDict("equationOfState").lookup<scalar>("dRho"),
        dict.subDict("equationOfState").lookup<scalar>("mine"),
        dict.subDict("equationOfState").lookup<scalar>("de")
    ),
    TTable_
    (
        dict.subDict("thermodynamics").lookup<fileName>("file"),
        dict.subDict("thermodynamics").lookup<word>("mod"),
        dict.subDict("thermodynamics").lookup<word>("rhoMod"),
        dict.subDict("thermodynamics").lookup<word>("eMod"),
        dict.subDict("thermodynamics").lookup<label>("nRho"),
        dict.subDict("thermodynamics").lookup<label>("ne"),
        dict.subDict("thermodynamics").lookup<scalar>("minRho"),
        dict.subDict("thermodynamics").lookup<scalar>("dRho"),
        dict.subDict("thermodynamics").lookup<scalar>("mine"),
        dict.subDict("thermodynamics").lookup<scalar>("de")
    )
{}

// ************************************************************************* //
