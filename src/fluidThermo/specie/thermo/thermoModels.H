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

#ifndef thermoModels_H
#define thermoModels_H

#include "constTransport.H"

#include "thermoModel.H"

#include "hConst.H"
#include "eConst.H"
#include "janafThermo.H"

#include "equationOfState.H"
#include "MGEquationOfState.H"

#include "CochranChan.H"
#include "Murnaghan.H"
#include "BirchMurnaghan2.H"
#include "BirchMurnaghan3.H"
#include "solidJWL.H"

#include "idealGas.H"
#include "stiffenedGas.H"
#include "Tait.H"
#include "vanderWaals.H"
#include "JWL.H"
#include "JWLC.H"
#include "LSZK.H"
#include "BWR.H"
#include "BKW.H"
#include "DoanNickel.H"

#include "tabulatedThermo.H"
#include "tabulatedThermoEOS.H"

#include "blastSpecie.H"

namespace Foam
{
// Solid materials
    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        CochranChan
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        Murnaghan
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BirchMurnaghan2
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BirchMurnaghan3
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        solidJWL
    );

// Fluids
    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        idealGas
    );
    defineThermoType
    (
        constTransport,
        janaf,
        MGEquationOfState,
        idealGas
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BKW
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        BWR
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        stiffenedGas
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        Tait
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        vanderWaals
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        JWL
    );

    defineThermoTypes
    (
        constTransport,
        equationOfState,
        JWLC
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        LSZK
    );

    defineThermoTypes
    (
        constTransport,
        MGEquationOfState,
        DoanNickel
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
        > constTransporttabulatedMGEquationOfStateDoanNickel;

    typedef constTransport<tabulatedThermoEOS<blastSpecie>>
        constTransporttabulatedMGEquationOfStatetabulated;
}

#endif
// ************************************************************************* //
