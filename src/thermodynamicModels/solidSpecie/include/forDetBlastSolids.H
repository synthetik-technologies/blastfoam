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

#ifndef forDetBlastSolids_H
#define forDetBlastSolids_H

#include "constIsoSolidTransport.H"
#include "constAnIsoSolidTransport.H"

#include "thermoModel.H"
#include "eConstBlastThermo.H"
#include "hConstBlastThermo.H"
#include "eTabulatedBlastThermo.H"
#include "hTabulatedBlastThermo.H"

#include "rhoConst.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forDetBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forDetrSolidEqns(uMu, rMu, uCp, rCp, uEos, Macro, Args...)            \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, rhoConst, specieBlast, Macro, Args);

#define forDetuSolidEqns(uMu, rMu, uCp, rCp, Macro, Args...)                  \
    forDetrSolidEqns(uMu, rMu, uCp, rCp, rhoConst, Macro, Args);

#define forDetrSolidThermos(uMu, rMu, uCp, Macro, Args...)                    \
    forDetuSolidEqns(uMu, rMu, uCp, eConstThermo, Macro, Args);               \
    forDetuSolidEqns(uMu, rMu, uCp, hConstThermo, Macro, Args);               \
    forDetuSolidEqns(uMu, rMu, uCp, eTabulatedThermo, Macro, Args);

#define forDetuSolidThermos(uMu, rMu, Macro, Args...)                         \
    forDetrSolidThermos(uMu, rMu, eConstThermo, Macro, Args);                 \
    forDetrSolidThermos(uMu, rMu, hConstThermo, Macro, Args);                 \
    forDetrSolidThermos(uMu, rMu, eTabulatedThermo, Macro, Args);

#define forDetrSolidTransports(uMu, Macro, Args...)                           \
    forDetuSolidThermos(uMu, constIsoSolidTransport, Macro, Args);            \
    forDetuSolidThermos(uMu, constAnIsoSolidTransport, Macro, Args);

#define forDetuSolidTransports(Macro, Args...)                                \
    forDetrSolidTransports(constIsoSolidTransport, Macro, Args);              \
    forDetrSolidTransports(constAnIsoSolidTransport, Macro, Args);

#define forDetSolids(Macro, Args...)                                         \
    forDetuSolidTransports(Macro, Args)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
// ************************************************************************* //
