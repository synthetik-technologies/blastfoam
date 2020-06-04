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

#include "solidThermoModel.H"
#include "basicSolidThermoModel.H"
#include "eThermoModel.H"

#include "constSolidIsoTransport.H"
#include "constSolidAnIsoTransport.H"

#include "thermoModel.H"
#include "eConst.H"
#include "hConst.H"

#include "equationOfState.H"
#include "rhoConst.H"

#include "blastSpecie.H"

#include "addToRunTimeSelectionTable.H"

namespace Foam
{
    addSolidThermo
    (
        constSolidIsoTransport,
        eConst,
        equationOfState,
        rhoConst
    );

    addSolidThermo
    (
        constSolidAnIsoTransport,
        eConst,
        equationOfState,
        rhoConst
    );

    addSolidThermo
    (
        constSolidIsoTransport,
        hConst,
        equationOfState,
        rhoConst
    );

    addSolidThermo
    (
        constSolidAnIsoTransport,
        hConst,
        equationOfState,
        rhoConst
    );


}
// ************************************************************************* //
