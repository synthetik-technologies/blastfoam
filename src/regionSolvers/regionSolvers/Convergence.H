/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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


Description
    Enumerations for iteration type and convergence

SourceFiles
    regionSolverList.C

\*---------------------------------------------------------------------------*/

#ifndef Convergence_H
#define Convergence_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// Iteration types
enum IterType
{
    // First iteration
    FIRST_ITER,

    // Not first or last iteration
    MID_ITER,

    // Final iteration
    FINAL_ITER
};


// Convergence types
enum Convergence
{
    // No convergence criteria has been specified
    UNKNOWN_CONVERGENCE,

    // Criteria not met
    NOT_CONVERGED,

    // Criteria met, but not 0
    CONVERGED,

    // Negligible error, assume no clean up steps are needed
    FULL_CONVERGENCE
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
