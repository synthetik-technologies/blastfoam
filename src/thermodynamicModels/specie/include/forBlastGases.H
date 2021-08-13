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

#ifndef forBlastGasses_H
#define forBlastGasses_H

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermoModel.H"

#include "hConstBlastThermo.H"
#include "eConstBlastThermo.H"
#include "janafBlastThermo.H"
#include "tabulatedBlastThermo.H"

#include "perfectGas.H"
#include "AbelNobel.H"

#include "idealGas.H"
#include "stiffenedGas.H"
#include "vanderWaals.H"
#include "LSZK.H"
#include "BKW.H"
#include "DoanNickel.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forGasEqns(Mu, Cp, Macro, Args...)                                  \
    forThermo(Mu, Cp, perfectGas, specieBlast, Macro, Args);                \
    forThermo(Mu, Cp, AbelNobel, specieBlast, Macro, Args);                 \
    forThermo(Mu, Cp, idealGas, specieBlast, Macro, Args);                  \
    forThermo(Mu, Cp, stiffenedGas, specieBlast, Macro, Args);              \
    forThermo(Mu, Cp, vanderWaals, specieBlast, Macro, Args);               \
    forThermo(Mu, Cp, LSZK, specieBlast, Macro, Args);                      \
    forThermo(Mu, Cp, BKW, specieBlast, Macro, Args);                       \
    forThermo(Mu, Cp, DoanNickel, specieBlast, Macro, Args)

#define forGasThermos(Mu, Macro, Args...)                                   \
    forGasEqns(Mu, eConstThermo, Macro, Args);                              \
    forGasEqns(Mu, janafThermo, Macro, Args);                               \
    forGasEqns(Mu, hConstThermo, Macro, Args)

#define forGasTransports(Macro, Args...)                                    \
    forGasThermos(constTransport, Macro, Args);                             \
    forGasThermos(sutherlandTransport, Macro, Args)

#define forGases(Macro, Args...)                                            \
    forGasTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCoeffGasEqns(Mu, Cp, Macro, Args...)                             \
    forThermo(Mu, Cp, perfectGas, specieBlast, Macro, Args);                \
    forThermo(Mu, Cp, AbelNobel, specieBlast, Macro, Args)

#define forCoeffGasThermos(Mu, Macro, Args...)                              \
    forCoeffGasEqns(Mu, hConstThermo, Macro, Args);                         \
    forCoeffGasEqns(Mu, janafThermo, Macro, Args);                          \
    forCoeffGasEqns(Mu, eConstThermo, Macro, Args)

#define forCoeffGasTransports(Macro, Args...)                               \
    forCoeffGasThermos(constTransport, Macro, Args);                        \
    forCoeffGasThermos(sutherlandTransport, Macro, Args)

#define forCoeffGases(Macro, Args...)                                       \
    forCoeffGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
