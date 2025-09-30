/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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

#ifndef makeMultivariateIntegrators_H
#define makeMultivariateIntegrators_H

#include "MultivariateIntegrator.H"
#include "MidPointMultivariateIntegrator.H"
#include "TrapezoidalMultivariateIntegrator.H"
#include "GaussianMultivariateIntegrator.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Define all integrators for a given type
#define makeMultivariateIntegratorTypes(Type, Adapt)                                  \
    makeMultivariateIntegratorType(MidPointMultivariateIntegrator, Type, Adapt)       \
    makeMultivariateIntegratorType(TrapezoidalMultivariateIntegrator, Type, Adapt)    \
    makeMultivariateIntegratorType(GaussianMultivariateIntegrator, Type, Adapt)




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Make the base integrator type along with the actual models
#define makeMultivariateIntegrators(Type, Adapt)                                      \
    makeMultivariateIntegrator(Type, Adapt)                                           \
    makeMultivariateIntegratorTypes(Type, Adapt);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Macro with a dummy first input for use with some macro based construction
#define makeMultivariateIntegratorsArg(Dummy, Type, Adapt)                            \
    makeMultivariateIntegrator(Type, Adapt)                                           \
    makeMultivariateIntegratorTypes(Type, Adapt);

#endif


// ************************************************************************* //
