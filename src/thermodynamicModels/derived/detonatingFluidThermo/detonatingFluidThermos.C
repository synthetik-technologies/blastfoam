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
        idealGas,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        idealGas,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        idealGas,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        idealGas,
        BKW
    );

//- Stiffened gas
    addDetonatingFluidThermos
    (
        constTransport,
        stiffenedGas,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        stiffenedGas,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        stiffenedGas,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        stiffenedGas,
        BKW
    );

//- Murnaghan
    addDetonatingFluidThermos
    (
        constTransport,
        Murnaghan,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        Murnaghan,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        Murnaghan,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        Murnaghan,
        BKW
    );

//- Birch-Murnaghan 2
    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan2,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan2,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan2,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan2,
        BKW
    );


//- Birch-Murnaghan 3
    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan3,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan3,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan3,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        BirchMurnaghan3,
        BKW
    );

//- Cochran Chan
    addDetonatingFluidThermos
    (
        constTransport,
        CochranChan,
        JWL
    );

    addDetonatingFluidThermos
    (
        constTransport,
        CochranChan,
        JWLC
    );

    addDetonatingFluidThermos
    (
        constTransport,
        CochranChan,
        LSZK
    );

    addDetonatingFluidThermos
    (
        constTransport,
        CochranChan,
        BKW
    );

// Blended JWL
    addDetonatingFluidThermos
    (
        constTransport,
        solidJWL,
        JWL
    );
}
// ************************************************************************* //
