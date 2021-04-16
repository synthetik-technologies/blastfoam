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

#ifndef fluidThermoModelTypes_H
#define fluidThermoModelTypes_H

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "thermoModel.H"

#include "hConst.H"
#include "eConst.H"
#include "janaf.H"

#include "equationOfState.H"
#include "MGEquationOfState.H"

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
#include "Tillotson.H"

#include "tabulatedThermo.H"
#include "tabulatedThermoEOS.H"

#include "blastSpecie.H"
#include "rBlastSpecie.H"

namespace Foam
{
// Solid materials
    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        CochranChan,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        Murnaghan,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BirchMurnaghan2,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        solidJWL,
        blastSpecie
    );

// Fluids
    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        idealGas,
        blastSpecie
    );
    defineThermoTypes
    (
        sutherlandTransport,
        MGEquationOfState,
        idealGas,
        blastSpecie
    );
    defineThermoType
    (
        constTransport,
        janaf,
        MGEquationOfState,
        idealGas,
        blastSpecie
    );
    defineThermoType
    (
        sutherlandTransport,
        janaf,
        MGEquationOfState,
        idealGas,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        perfectGas,
        blastSpecie
    );
    defineThermoTypes
    (
        sutherlandTransport,
        equationOfState,
        perfectGas,
        blastSpecie
    );
    defineThermoType
    (
        constTransport,
        janaf,
        equationOfState,
        perfectGas,
        blastSpecie
    );
    defineThermoType
    (
        sutherlandTransport,
        janaf,
        equationOfState,
        perfectGas,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        AbelNobel,
        blastSpecie
    );
    defineThermoTypes
    (
        sutherlandTransport,
        equationOfState,
        AbelNobel,
        blastSpecie
    );
    defineThermoTypes
    (
        propellantTransport,
        equationOfState,
        AbelNobel,
        blastSpecie
    );
    defineThermoType
    (
        constTransport,
        janaf,
        equationOfState,
        AbelNobel,
        blastSpecie
    );
    defineThermoType
    (
        sutherlandTransport,
        janaf,
        equationOfState,
        AbelNobel,
        blastSpecie
    );
    defineThermoType
    (
        propellantTransport,
        janaf,
        equationOfState,
        AbelNobel,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BKW,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BWR,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        Tait,
        blastSpecie
    );

    typedef
        constTransport
        <
            thermoModel
            <
                eConst
                <
                    Tillotson<blastSpecie>
                >
            >
        > constTransporteConstTillotsonblastSpecie;

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        vanderWaals,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        JWL,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        JWLC,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        LSZK,
        blastSpecie
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        DoanNickel,
        blastSpecie
    );

    typedef
        constTransport
        <
            tabulatedThermo
            <
                MGEquationOfState
                <
                    DoanNickel<blastSpecie>
                >
            >
        > constTransporttabulatedMGEquationOfStateDoanNickelblastSpecie;

    typedef constTransport<tabulatedThermoEOS<blastSpecie>>
        constTransporttabulatedMGEquationOfStatetabulatedblastSpecie;
}

#endif
// ************************************************************************* //
