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

#ifndef forDetBlastGasses_H
#define forDetBlastGasses_H

#include "constTransport.H"

#include "thermoModel.H"

#include "eConstBlastThermo.H"
#include "ePolynomialBlastThermo.H"

#include "perfectGas.H"
#include "AbelNobel.H"
#include "JWLC.H"

#include "vanderWaals.H"
#include "JWL.H"
#include "LSZK.H"
#include "BKW.H"

#include "CochranChan.H"
#include "Murnaghan.H"
#include "BirchMurnaghan2.H"
#include "BirchMurnaghan3.H"
#include "solidJWL.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forDetBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forDetrGasEqns(uMu, rMu, uCp, rCp, uEos, Macro, Args...)            \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, JWLC, specieBlast, Macro, Args); \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, JWL, specieBlast, Macro, Args);  \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, LSZK, specieBlast, Macro, Args); \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, BKW, specieBlast, Macro, Args);


#define forDetuGasEqns(uMu, rMu, uCp, rCp, Macro, Args...)                  \
    forDetrGasEqns(uMu, rMu, uCp, rCp, CochranChan, Macro, Args);           \
    forDetrGasEqns(uMu, rMu, uCp, rCp, Murnaghan, Macro, Args);             \
    forDetrGasEqns(uMu, rMu, uCp, rCp, BirchMurnaghan2, Macro, Args);       \
    forDetrGasEqns(uMu, rMu, uCp, rCp, BirchMurnaghan3, Macro, Args);       \
    forDetrGasEqns(uMu, rMu, uCp, rCp, solidJWL, Macro, Args);

#define forDetrGasThermos(uMu, rMu, uCp, Macro, Args...)                    \
    forDetuGasEqns(uMu, rMu, uCp, eConstThermo, Macro, Args);               \
    forDetuGasEqns(uMu, rMu, uCp, ePolynomialThermo, Macro, Args);

#define forDetuGasThermos(uMu, rMu, Macro, Args...)                         \
    forDetrGasThermos(uMu, rMu, eConstThermo, Macro, Args);                 \
    forDetrGasThermos(uMu, rMu, ePolynomialThermo, Macro, Args);

#define forDetrGasTransports(uMu, Macro, Args...)                           \
    forDetuGasThermos(uMu, constTransport, Macro, Args);

#define forDetuGasTransports(Macro, Args...)                                \
    forDetrGasTransports(constTransport, Macro, Args);

#define forDetGases(Macro, Args...)                                         \
    forDetuGasTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
