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

#include "phaseFluidBlastThermo.H"
#include "cavitatingFluidBlastThermo.H"
#include "blendedBlastThermo.H"
#include "forDetBlastGases.H"
#include "makeDetBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "constTransport.H"

#include "thermoModel.H"

#include "eConstBlastThermo.H"
#include "ePolynomialBlastThermo.H"

#include "idealGas.H"
#include "perfectGas.H"

#include "linearTillotson.H"

#include "specieBlast.H"
#include "rspecieBlast.H"

#include "forDetBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define forCavvFluidEqns(uMu, rMu, uCp, rCp, uEos, Macro, Args...)            \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, idealGas, specieBlast, Macro, Args); \
    forDetThermo(uMu, rMu, uCp, rCp, uEos, perfectGas, specieBlast, Macro, Args);


#define forCavlFluidEqns(uMu, rMu, uCp, rCp, Macro, Args...)                  \
    forCavvFluidEqns(uMu, rMu, uCp, rCp, linearTillotson, Macro, Args);

#define forCavvrFluidThermos(uMu, rMu, uCp, Macro, Args...)                    \
    forCavlFluidEqns(uMu, rMu, uCp, eConstThermo, Macro, Args);

#define forCavlFluidThermos(uMu, rMu, Macro, Args...)                         \
    forCavvrFluidThermos(uMu, rMu, eConstThermo, Macro, Args);

#define forCavvFluidTransports(uMu, Macro, Args...)                           \
    forCavlFluidThermos(uMu, constTransport, Macro, Args);

#define forCavlGasTransports(Macro, Args...)                                \
    forCavvFluidTransports(constTransport, Macro, Args);

#define forCavFluids(Macro, Args...)                                         \
    forCavlGasTransports(Macro, Args)

namespace Foam
{
    forCavFluids
    (
        makeDetThermo,
        phaseFluidBlastThermo,
        cavitatingFluidBlastThermo,
        blendedBlastThermo
    );
}
// ************************************************************************* //
