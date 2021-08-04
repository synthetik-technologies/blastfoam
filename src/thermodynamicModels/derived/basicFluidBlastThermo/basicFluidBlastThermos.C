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
#include "basicFluidBlastThermo.H"
#include "eBlastThermo.H"
#include "forBlastGases.H"
#include "makeBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
//     typedefThermo(constTransport, eConst, perfectGas, blastSpecie);
//     typedef
//         basicFluidBlastThermo
//         <
//             eBlastThermo
//             <
//                 phaseFluidBlastThermo, constTransporteConstperfectGasblastSpecie
//             >
//         > basicFluidBlastThermoconstTransporteConstperfectGasblastSpecie;
    forGases
    (
        makeThermo,
        phaseFluidBlastThermo,
        basicFluidBlastThermo,
        eBlastThermo
    );
// BaseThermo = phaseFluidBlastThermo, etc
// DerivThermo = basicFluidBlastThermo, etc
// CThermo = eThermo or blendedThermo
// ThermoPhys = Transport##Thermo##Eos##Specie
// Solids
//     addFluidThermo
//     (
//         phaseFluidBlastThermo,
//         basicFluidBlastThermo,
//         constTransport,
//         eConst,
//         CochranChan
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         Murnaghan
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         BirchMurnaghan2
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         BirchMurnaghan3
//     );
//
// // Fluids
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         idealGas
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         perfectGas
//     );
//     addFluidThermos
//     (
//         sutherlandTransport,
//         perfectGas
//     );
//     addFluidThermo
//     (
//         sutherlandTransport,
//         janaf,
//         perfectGas
//     );
//     addFluidThermo
//     (
//         constTransport,
//         janaf,
//         perfectGas
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         stiffenedGas
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         Tait
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         linearTillotson
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         Tillotson
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         AbelNobel
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         BKW
//     );
//
//     addFluidThermos
//     (
//         constTransport,
//         BWR
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         vanderWaals
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         JWL
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         JWLC
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         LSZK
//     );
//
//     addFluidThermo
//     (
//         constTransport,
//         eConst,
//         DoanNickel
//     );
//     addFluidThermo
//     (
//         constTransport,
//         tabulatedThermo,
//         DoanNickel
//     );
//
//
//     // Tabulated equation of state and thermo model
//     typedef basicFluidBlastThermo
//         <
//             eBlastThermo
//             <
//                 fluidBlastThermo,
//                 constTransporttabulatedMGEquationOfStatetabulatedblastSpecie
//             >
//         >
//         basicFluidThermoconstTransporttabulatedMGEquationOfStatetabulated;
//
//     defineTemplateTypeNameAndDebugWithName
//     (
//         basicFluidThermoconstTransporttabulatedMGEquationOfStatetabulated,
//         (constTransporttabulatedMGEquationOfStatetabulatedblastSpecie::typeName()).c_str(),
//         0
//     );
//     addToRunTimeSelectionTable
//     (
//         fluidBlastThermo,
//         basicFluidThermoconstTransporttabulatedMGEquationOfStatetabulated,
//         basic
//     );
}
// ************************************************************************* //
