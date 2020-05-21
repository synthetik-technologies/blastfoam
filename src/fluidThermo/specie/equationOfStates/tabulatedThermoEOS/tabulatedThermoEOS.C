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
        dict.subDict("equationOfState").lookupType<fileName>("file"),
        dict.subDict("equationOfState").lookup("mod"),
        dict.subDict("equationOfState").lookup("rhoMod"),
        dict.subDict("equationOfState").lookup("eMod"),
        dict.subDict("equationOfState").lookupType<label>("nRho"),
        dict.subDict("equationOfState").lookupType<label>("ne"),
        dict.subDict("equationOfState").lookupType<scalar>("minRho"),
        dict.subDict("equationOfState").lookupType<scalar>("dRho"),
        dict.subDict("equationOfState").lookupType<scalar>("mine"),
        dict.subDict("equationOfState").lookupType<scalar>("de")
    ),
    TTable_
    (
        dict.subDict("thermodynamics").lookupType<fileName>("file"),
        dict.subDict("thermodynamics").lookup("mod"),
        dict.subDict("thermodynamics").lookup("rhoMod"),
        dict.subDict("thermodynamics").lookup("eMod"),
        dict.subDict("thermodynamics").lookupType<label>("nRho"),
        dict.subDict("thermodynamics").lookupType<label>("ne"),
        dict.subDict("thermodynamics").lookupType<scalar>("minRho"),
        dict.subDict("thermodynamics").lookupType<scalar>("dRho"),
        dict.subDict("thermodynamics").lookupType<scalar>("mine"),
        dict.subDict("thermodynamics").lookupType<scalar>("de")
    )
{}

// ************************************************************************* //
