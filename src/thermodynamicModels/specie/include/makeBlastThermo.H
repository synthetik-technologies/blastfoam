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

#ifndef makeBlastThermo_H
#define makeBlastThermo_H

#include "basicBlastThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineThermo(BaseThermo, DerivThermo, CThermo, ThermoPhys)          \
    typedef                                                                 \
        DerivThermo                                                         \
        <                                                                   \
            CThermo                                                         \
            <                                                               \
                BaseThermo, ThermoPhys                                      \
            >                                                               \
        > DerivThermo##ThermoPhys;                                          \
                                                                            \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        DerivThermo##ThermoPhys,                                            \
        (                                                                   \
            word(BaseThermo::typeName_()) + "<"                             \
          + word(DerivThermo::typeName_()) + "<"                            \
          + ThermoPhys::typeName() + ">>"                                   \
        ).c_str(),                                                          \
        0                                                                   \
    );

#define makeThermo(BaseThermo, DerivThermo, CThermo, ThermoPhys)            \
                                                                            \
    defineThermo(BaseThermo, DerivThermo, CThermo, ThermoPhys);             \
                                                                            \
    addThermo(BaseThermo, DerivThermo, CThermo, ThermoPhys);

#define addThermo(BaseThermo, DerivThermo, CThermo, ThermoPhys)             \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        BaseThermo,                                                         \
        DerivThermo##ThermoPhys,                                            \
        dictionary                                                          \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
