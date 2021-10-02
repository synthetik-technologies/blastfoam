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

#ifndef forBlastSolidFluids_H
#define forBlastSolidFluids_H

#include "constTransport.H"

#include "thermoModel.H"

#include "eConstBlastThermo.H"
#include "ePolynomialBlastThermo.H"

#include "CochranChan.H"
#include "Murnaghan.H"
#include "BirchMurnaghan2.H"
#include "BirchMurnaghan3.H"
#include "solidJWL.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forSolidFluidEqns(Mu, Cp, Macro, Args...)                           \
    forThermo(Mu, Cp, CochranChan, specieBlast, Macro, Args);               \
    forThermo(Mu, Cp, Murnaghan, specieBlast, Macro, Args);                 \
    forThermo(Mu, Cp, BirchMurnaghan2, specieBlast, Macro, Args);           \
    forThermo(Mu, Cp, BirchMurnaghan3, specieBlast, Macro, Args);           \
    forThermo(Mu, Cp, solidJWL, specieBlast, Macro, Args)

#define forSolidFluidThermos(Mu, Macro, Args...)                            \
    forSolidFluidEqns(Mu, eConstThermo, Macro, Args);                       \
    forSolidFluidEqns(Mu, ePolynomialThermo, Macro, Args);

#define forSolidFluidTransports(Macro, Args...)                             \
    forSolidFluidThermos(constTransport, Macro, Args);

#define forSolidFluids(Macro, Args...)                                      \
    forSolidFluidTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
