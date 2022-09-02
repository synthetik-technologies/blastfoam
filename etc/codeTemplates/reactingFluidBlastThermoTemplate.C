/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "forBlastThermo.H"
#include "makeBlastThermo.H"

#include "${specie}.H"

#include "thermoModel.H"

// EoS
#include "${equationOfState}.H"

// Thermo
#include "${thermo}BlastThermo.H"

// Transport
#include "${transport}Transport.H"

// reacting
#include "fluidBlastThermo.H"
#include "reactingFluidBlastThermo.H"
#include "mixtureBlastThermo.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ${SHA1sum}
    //
    // Unique function name that can be checked if the correct library version
    // has been loaded
    void ${typeName}_${SHA1sum}(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define ThermoPhysics                                                          \
    ${transport}Transport${thermo}Thermo${equationOfState}${specie}

namespace Foam
{
    forThermo
    (
        ${transport}Transport,
        ${thermo}Thermo,
        ${equationOfState},
        ${specie},
        makeThermo,
        fluidBlastThermo,
        reactingFluidBlastThermo,
        mixtureBlastThermo
    );

    // Define typeName and debug for multicomponentFLuidBlastThermo
    forThermo
    (
        ${transport}Transport,
        ${thermo}Thermo,
        ${equationOfState},
        ${specie},
        defineThermo,
        fluidBlastThermo,
        multicomponentFluidBlastThermo,
        mixtureBlastThermo
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "noBlastChemistrySolver.H"
#include "EulerImplicitBlastChemistrySolver.H"
#include "odeBlastChemistrySolver.H"

#include "standardBlastChemistryModel.H"

#include "makeBlastChemistrySolver.H"

namespace Foam
{
    defineChemistrySolvers(fluidBlastThermo, ThermoPhysics);

    makeChemistrySolvers(noBlastChemistrySolver, fluidBlastThermo, ThermoPhysics);
    makeChemistrySolvers(EulerImplicitBlastChemistrySolver, fluidBlastThermo, ThermoPhysics);
    makeChemistrySolvers(odeBlastChemistrySolver, fluidBlastThermo, ThermoPhysics);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "FallOffReactionRate.H"

#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"

#include "LangmuirHinshelwoodReactionRate.H"

namespace Foam
{
    defineReaction(nullArg, ThermoPhysics);

    makeIRNReactions(ArrheniusReactionRate, ThermoPhysics);
    makeIRNReactions(LandauTellerReactionRate, ThermoPhysics);
    makeIRNReactions(thirdBodyArrheniusReactionRate, ThermoPhysics);

    makeIRReactions(JanevReactionRate, ThermoPhysics);
    makeIRReactions(powerSeriesReactionRate, ThermoPhysics);
}

// ************************************************************************* //
