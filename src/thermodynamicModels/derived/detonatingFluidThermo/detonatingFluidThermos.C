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

#include "fluidThermoModel.H"
#include "blendedThermoModel.H"
#include "detonatingFluidThermo.H"
#include "thermoModels.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
//- Ideal gas
    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        constTransport,
        MGEquationOfState,
        JWL
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         MGEquationOfState,
//         idealGas,
//         constTransport,
//         equationOfState,
//         JWLC
//     );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        constTransport,
        MGEquationOfState,
        LSZK
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         MGEquationOfState,
//         idealGas,
//         constTransport,
//         equationOfState,
//         BKW
//     );

//- Stiffened gas
    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        constTransport,
        MGEquationOfState,
        JWL
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         MGEquationOfState,
//         stiffenedGas,
//         constTransport,
//         equationOfState,
//         JWLC
//     );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        constTransport,
        MGEquationOfState,
        LSZK
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         MGEquationOfState,
//         stiffenedGas,
//         constTransport,
//         equationOfState,
//         BKW
//     );

//- Murnaghan
    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan,
        constTransport,
        MGEquationOfState,
        JWL
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         Murnaghan,
//         constTransport,
//         equationOfState,
//         JWLC
//     );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan,
        constTransport,
        MGEquationOfState,
        LSZK
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         Murnaghan,
//         constTransport,
//         equationOfState,
//         BKW
//     );

// //- Birch-Murnaghan 2
//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan2,
//         constTransport,
//         MGEquationOfState,
//         JWL
//     );
//
//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan2,
//         constTransport,
//         equationOfState,
//         JWLC
//     );
//
//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan2,
//         constTransport,
//         MGEquationOfState,
//         LSZK
//     );
//
//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan2,
//         constTransport,
//         equationOfState,
//         BKW
//     );


//- Birch-Murnaghan 3
    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        constTransport,
        MGEquationOfState,
        JWL
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan3,
//         constTransport,
//         equationOfState,
//         JWLC
//     );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        constTransport,
        MGEquationOfState,
        LSZK
    );

//     addDetonatingFluidThermos
//     (
//         constTransport,
//         equationOfState,
//         BirchMurnaghan3,
//         constTransport,
//         equationOfState,
//         BKW
//     );

// Blended JWL
    addDetonatingFluidThermo
    (
        constTransport,
        hConst,
        MGEquationOfState,
        solidJWL,
        constTransport,
        eConst,
        MGEquationOfState,
        JWL
    );
}
// ************************************************************************* //
