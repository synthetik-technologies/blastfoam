/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "makeReaction.H"

#include "ArrheniusReactionRate.H"
#include "infiniteReactionRate.H"
#include "LandauTellerReactionRate.H"
#include "thirdBodyArrheniusReactionRate.H"

#include "ChemicallyActivatedReactionRate.H"
#include "JanevReactionRate.H"
#include "powerSeriesReactionRate.H"

#include "FallOffReactionRate.H"
#include "LindemannFallOffFunction.H"
#include "SRIFallOffFunction.H"
#include "TroeFallOffFunction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define makeReactionsFromTypes(Transport, Thermo, EOS)                      \
    makeReactions                                                           \
    (                                                                       \
        Transport##Thermo##EOS##blastSpecie,                                \
        Reaction<Transport##Thermo##EOS##blastSpecie>                       \
    )


#define makeReactions(Thermo, Reaction)                                        \
                                                                               \
    defineTemplateTypeNameAndDebug(Reaction, 0);                               \
    defineTemplateRunTimeSelectionTable(Reaction, dictionary);                 \
                                                                               \
    makeIRNReactions(Thermo, ArrheniusReactionRate)                            \
    makeIRNReactions(Thermo, infiniteReactionRate)                             \
    makeIRNReactions(Thermo, LandauTellerReactionRate)                         \
    makeIRNReactions(Thermo, thirdBodyArrheniusReactionRate)                   \
                                                                               \
    makeIRReactions(Thermo, JanevReactionRate)                                 \
    makeIRReactions(Thermo, powerSeriesReactionRate)                           \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       ArrheniusReactionRate,                                                  \
       LindemannFallOffFunction                                                \
    )                                                                          \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       ArrheniusReactionRate,                                                  \
       TroeFallOffFunction                                                     \
    )                                                                          \
                                                                               \
    makePressureDependentReactions                                             \
    (                                                                          \
       Thermo,                                                                 \
       ArrheniusReactionRate,                                                  \
       SRIFallOffFunction                                                      \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // sensible enthalpy based reactions
    makeReactionsFromTypes(constTransport, eConst, idealGas)
    makeReactionsFromTypes(constTransport, hConst, idealGas)
    makeReactionsFromTypes(constTransport, janaf, idealGas)
    makeReactionsFromTypes(sutherlandTransport, eConst, idealGas)
    makeReactionsFromTypes(sutherlandTransport, hConst, idealGas)
    makeReactionsFromTypes(sutherlandTransport, janaf, idealGas)

    makeReactionsFromTypes(constTransport, eConst, perfectGas)
    makeReactionsFromTypes(constTransport, hConst, perfectGas)
    makeReactionsFromTypes(constTransport, janaf, perfectGas)
    makeReactionsFromTypes(sutherlandTransport, eConst, perfectGas)
    makeReactionsFromTypes(sutherlandTransport, hConst, perfectGas)
    makeReactionsFromTypes(sutherlandTransport, janaf, perfectGas)
}

// ************************************************************************* //
