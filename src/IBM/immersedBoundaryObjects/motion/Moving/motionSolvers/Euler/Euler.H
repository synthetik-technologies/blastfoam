/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2019 OpenFOAM Foundation
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

Class
    Foam::immersedBoundaryObjectMotionSolvers::Euler

Description
    First order Euler time integration

    Example specification in dynamicMeshDict:
    \verbatim
    solver
    {
        type    Euler;
    }
    \endverbatim

See also
    Foam::immersedBoundaryObjectMotionSolvers::Euler

SourceFiles
    Euler.C

\*---------------------------------------------------------------------------*/

#ifndef EulerMotionSolver_H
#define EulerMotionSolver_H

#include "immersedBoundaryObjectMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace immersedBoundaryObjectMotionSolvers
{

/*---------------------------------------------------------------------------*\
                           Class Euler Declaration
\*---------------------------------------------------------------------------*/

class Euler
:
    public immersedBoundaryObjectMotionSolver
{
public:

    //- Runtime type information
    TypeName("Euler");


    // Constructors

        //- Construct from a dictionary and the body
        Euler
        (
            const dictionary& dict,
            immersedBoundaryObject& body,
            movingImmersedBoundaryObjectState& state,
            movingImmersedBoundaryObjectState& state0
        );


    //- Destructor
    virtual ~Euler();


    // Member Functions

        //- Drag coefficient
        virtual void solve
        (
            bool firstIter,
            const vector& fGlobal,
            const vector& tauGlobal,
            scalar deltaT,
            scalar deltaT0
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace immersedBoundaryObjectMotionSolvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
