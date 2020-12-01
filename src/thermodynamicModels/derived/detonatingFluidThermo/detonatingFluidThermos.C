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
#include "fluidThermoModelTypes.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
//- Ideal gas
    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        MGEquationOfState,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        equationOfState,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        MGEquationOfState,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        equationOfState,
        BKW
    );

//- Stiffened gas
    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        MGEquationOfState,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        equationOfState,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        MGEquationOfState,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        equationOfState,
        BKW
    );

//- Murnaghan
    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan,
        MGEquationOfState,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan,
        equationOfState,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan,
        MGEquationOfState,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        Murnaghan,
        equationOfState,
        BKW
    );

//- Birch-Murnaghan 2
    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan2,
        MGEquationOfState,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan2,
        equationOfState,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan2,
        MGEquationOfState,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan2,
        equationOfState,
        BKW
    );


//- Birch-Murnaghan 3
    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        MGEquationOfState,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        equationOfState,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        MGEquationOfState,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        equationOfState,
        BKW
    );

//- Cochran Chan
    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        CochranChan,
        MGEquationOfState,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        CochranChan,
        equationOfState,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        CochranChan,
        MGEquationOfState,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        CochranChan,
        equationOfState,
        BKW
    );

// Blended JWL
    addDetonatingFluidThermos
    (
        constTransport,
        MGEquationOfState,
        solidJWL,
        MGEquationOfState,
        JWL
    );
}
// ************************************************************************* //
