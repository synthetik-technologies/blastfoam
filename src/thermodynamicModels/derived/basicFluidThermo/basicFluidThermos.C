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

#include "basicFluidThermo.H"
#include "eThermoModel.H"
#include "thermoModels.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
// Solids
    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        CochranChan
    );

    addFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan
    );

//     addFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan2
//     );

    addFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3
    );

// Fluids
    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas
    );
    addFluidThermo
    (
        constTransport,
        janaf,
        MGEquationOfState,
        idealGas
    );

    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas
    );

    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        Tait
    );

//     addFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BKW
//     );

    addFluidThermos
    (
        constTransport,
        equationOfState,
        BWR
    );

    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        vanderWaals
    );

    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        JWL
    );

//     addFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         JWLC
//     );

    addFluidThermos
    (
        constTransport,
        MGEquationOfState,
        LSZK
    );

    addFluidThermo
    (
        constTransport,
        eConst,
        MGEquationOfState,
        DoanNickel
    );

    // Tabulated equation of state and thermo model
    typedef basicFluidThermo
        <
            eThermoModel
            <
                fluidThermoModel,
                constTransporttabulatedMGEquationOfStateDoanNickel
            >
        >
        basicFluidThermoconstTransporttabulatedMGEquationOfStateDoanNickel;

    defineTemplateTypeNameAndDebugWithName
    (
        basicFluidThermoconstTransporttabulatedMGEquationOfStateDoanNickel,
        (constTransporttabulatedMGEquationOfStateDoanNickel::typeName()).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        fluidThermoModel,
        basicFluidThermoconstTransporttabulatedMGEquationOfStateDoanNickel,
        basic
    );

    // Tabulated equation of state and thermo model
    typedef basicFluidThermo
        <
            eThermoModel
            <
                fluidThermoModel,
                constTransporttabulatedMGEquationOfStatetabulated
            >
        >
        basicFluidThermoconstTransporttabulatedMGEquationOfStatetabulated;

    defineTemplateTypeNameAndDebugWithName
    (
        basicFluidThermoconstTransporttabulatedMGEquationOfStatetabulated,
        (constTransporttabulatedMGEquationOfStatetabulated::typeName()).c_str(),
        0
    );
    addToRunTimeSelectionTable
    (
        fluidThermoModel,
        basicFluidThermoconstTransporttabulatedMGEquationOfStatetabulated,
        basic
    );
}
// ************************************************************************* //
