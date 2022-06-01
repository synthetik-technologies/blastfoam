/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#ifndef makeBlastChemistrySolver_H
#define makeBlastChemistrySolver_H

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolver(Model, Thermo, ThermoPhysics)                    \
                                                                               \
    typedef Model<Thermo, ThermoPhysics>                                       \
        Model##Thermo##ThermoPhysics;                                          \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Model##Thermo##ThermoPhysics,                                          \
        (                                                                      \
            word(Model##Thermo##ThermoPhysics::typeName_())                    \
          + "<" + ThermoPhysics::typeName() + ">"                              \
        ).c_str(),                                                             \
        0                                                                      \
    )


#define makeChemistrySolver(Solver, Thermo, Model, ThermoPhysics)              \
                                                                               \
    typedef Solver<Model<Thermo, ThermoPhysics>>                               \
        Solver##Model##Thermo##ThermoPhysics;                                  \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        Solver##Model##Thermo##ThermoPhysics,                                  \
        (                                                                      \
            word(Solver##Model##Thermo##ThermoPhysics::typeName_())            \
          + "<" + word(Model<Thermo, ThermoPhysics>::typeName_())              \
          + "<" + ThermoPhysics::typeName() + ">>"                             \
        ).c_str(),                                                             \
        0                                                                      \
    );                                                                         \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        basicBlastChemistryModel,                                              \
        Solver##Model##Thermo##ThermoPhysics,                                  \
        thermo                                                                 \
    )


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineChemistrySolvers(Thermo, ThermoPhysics)                          \
    defineChemistrySolver                                                      \
    (                                                                          \
        standardBlastChemistryModel,                                           \
        Thermo,                                                                \
        ThermoPhysics                                                          \
    );

#define makeChemistrySolvers(Solver, Thermo, ThermoPhysics)                    \
    makeChemistrySolver                                                        \
    (                                                                          \
        Solver,                                                                \
        Thermo,                                                                \
        standardBlastChemistryModel,                                           \
        ThermoPhysics                                                          \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
