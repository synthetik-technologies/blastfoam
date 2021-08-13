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

#ifndef forBlastSolids_H
#define forBlastSolids_H

#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"

#include "thermoModel.H"
#include "eConstBlastThermo.H"
#include "hConstBlastThermo.H"
#include "tabulatedBlastThermo.H"

#include "rhoConst.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forSolidEqns(Mu, Cp, Macro, Args...)                                \
    forThermo(Mu, Cp, rhoConst, specieBlast, Macro, Args);

#define forSolidThermos(Mu, Macro, Args...)                                 \
    forSolidEqns(Mu, eConstThermo, Macro, Args);                            \
    forSolidEqns(Mu, tabulatedThermo, Macro, Args);                         \
    forSolidEqns(Mu, hConstThermo, Macro, Args)

#define forSolidTransports(Macro, Args...)                                  \
    forSolidThermos(constIsoSolidTransport, Macro, Args);                   \
    forSolidThermos(constAnIsoSolidTransport, Macro, Args)

#define forSolids(Macro, Args...)                                           \
    forSolidTransports(Macro, Args)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
