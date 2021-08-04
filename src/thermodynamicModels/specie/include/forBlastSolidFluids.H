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

#ifndef fluidThermoModelTypes_H
#define fluidThermoModelTypes_H

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermoModel.H"

#include "hConst.H"
#include "eConst.H"
#include "janaf.H"

#include "equationOfState.H"

#include "CochranChan.H"
#include "Murnaghan.H"
#include "BirchMurnaghan2.H"
#include "BirchMurnaghan3.H"
#include "solidJWL.H"

#include "idealGas.H"
#include "perfectGas.H"
#include "stiffenedGas.H"
#include "Tait.H"
#include "AbelNobel.H"
#include "vanderWaals.H"
#include "JWL.H"
#include "JWLC.H"
#include "LSZK.H"
#include "BWR.H"
#include "BKW.H"
#include "DoanNickel.H"
#include "linearTillotson.H"
#include "Tillotson.H"

#include "tabulatedThermo.H"
#include "tabulatedThermoEOS.H"

#include "blastSpecie.H"
#include "rBlastSpecie.H"

#include "forBlastThermo.H"

namespace Foam
{
// Solid materials
    defineThermoType
    (
        constTransport,
        eConst,
        CochranChan,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        Murnaghan,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        BirchMurnaghan2,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        BirchMurnaghan3,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        solidJWL,
        blastSpecie
    );

// Fluids
    defineThermoType
    (
        constTransport,
        eConst,
        idealGas,
        blastSpecie
    );
    defineThermoType
    (
        sutherlandTransport,
        eConst,
        idealGas,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        perfectGas,
        blastSpecie
    );
    defineThermoTypes
    (
        sutherlandTransport,
        perfectGas,
        blastSpecie
    );
    defineThermoType
    (
        constTransport,
        janaf,
        perfectGas,
        blastSpecie
    );
    defineThermoType
    (
        sutherlandTransport,
        janaf,
        perfectGas,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        AbelNobel,
        blastSpecie
    );
    defineThermoTypes
    (
        sutherlandTransport,
        AbelNobel,
        blastSpecie
    );
    defineThermoType
    (
        constTransport,
        janaf,
        AbelNobel,
        blastSpecie
    );
    defineThermoType
    (
        sutherlandTransport,
        janaf,
        AbelNobel,
        blastSpecie
    );


    defineThermoTypes
    (
        constTransport,
        BKW,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        BWR,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        stiffenedGas,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        Tait,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        Tillotson,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        linearTillotson,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        vanderWaals,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        JWL,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        JWLC,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        LSZK,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        eConst,
        DoanNickel,
        blastSpecie
    );

    defineThermoType
    (
        constTransport,
        tabulatedThermo,
        DoanNickel,
        blastSpecie
    );

    typedef constTransport<tabulatedThermoEOS<blastSpecie>>
        constTransporttabulatedMGEquationOfStatetabulatedblastSpecie;
}

#endif
// ************************************************************************* //
