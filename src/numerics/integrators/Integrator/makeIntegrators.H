/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#ifndef makeIntegrators_H
#define makeIntegrators_H

#include "Integrator.H"
#include "MidPointIntegrator.H"
#include "TrapezoidalIntegrator.H"
#include "Simpson13Integrator.H"
#include "Simpson38Integrator.H"
#include "BooleIntegrator.H"
#include "GaussianIntegrator.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define all integrators for a given type
#define makeIntegratorTypes(Type, Adapt)                                       \
    makeIntegratorType(MidPointIntegrator, Type, Adapt)                        \
    makeIntegratorType(TrapezoidalIntegrator, Type, Adapt)                     \
    makeIntegratorType(Simpson13Integrator, Type, Adapt)                       \
    makeIntegratorType(Simpson38Integrator, Type, Adapt)                       \
    makeIntegratorType(GaussianIntegrator, Type, Adapt)                        \
    makeIntegratorType(BooleIntegrator, Type, Adapt)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Make the base integrator type along with the actual models
#define makeIntegrators(Type, Adapt)                                           \
    makeIntegrator(Type, Adapt)                                                \
    makeIntegratorTypes(Type, Adapt);

// Make the base integrator type along with the actual models
#define makeAdaptIntegrators(Adapt, Type)                                           \
    makeIntegrator(Type, Adapt)                                                \
    makeIntegratorTypes(Type, Adapt);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
