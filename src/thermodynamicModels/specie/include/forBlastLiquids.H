/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2021-08-02 Jeff Heylmun:    Added blastFoam thermo classes
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

#ifndef forBlastLiquids_H
#define forBlastLiquids_H

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermoModel.H"

#include "hConstBlastThermo.H"
#include "eConstBlastThermo.H"
#include "janafBlastThermo.H"

#include "Tait.H"
#include "BWR.H"
#include "linearTillotson.H"
#include "Tillotson.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forLiquidEqns(Mu, Cp, Macro, Args...)                               \
    forThermo(Mu, Cp, Tait, specieBlast, Macro, Args);                      \
    forThermo(Mu, Cp, BWR, specieBlast, Macro, Args);                       \
    forThermo(Mu, Cp, linearTillotson, specieBlast, Macro, Args);                      \
    forThermo(Mu, Cp, Tillotson, specieBlast, Macro, Args);                 \

#define forLiquidThermos(Mu, Macro, Args...)                                \
    forLiquidEqns(Mu, eConstThermo, Macro, Args);                           \
    forLiquidEqns(Mu, janafThermo, Macro, Args);                            \
    forLiquidEqns(Mu, hConstThermo, Macro, Args)

#define forLiquidTransports(Macro, Args...)                                 \
    forLiquidThermos(constTransport, Macro, Args);

#define forLiquids(Macro, Args...)                                          \
    forLiquidTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
