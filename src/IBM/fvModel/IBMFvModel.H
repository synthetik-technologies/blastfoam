/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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
    Foam::fv::IBMForce

Description
    Adds immersed boundary forcing

Usage
    Example usage:
    \verbatim
    buoyancyForce1
    {
        type        IBM;
    }
    \endverbatim

SourceFiles
    IBMForce.C

\*---------------------------------------------------------------------------*/

#ifndef IBMForceModel_H
#define IBMForceModel_H

#include "fvModel.H"
#include "fvConstraint.H"
#include "immersedBoundaryObjectListSolver.H"
#include "hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                        Class IBMForce Declaration
\*---------------------------------------------------------------------------*/

class IBMForceModel
:
    public fvModel
{
    // Private Data

        //- Immersed objects
        immersedBoundaryObjectListSolver& ibm_;


public:

    //- Runtime type information
    TypeName("IBM");


    // Constructors

        //- Construct from explicit source name and mesh
        IBMForceModel
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        IBMForceModel(const IBMForceModel&) = delete;


    // Member Functions

        //- Correct the fvModel
        //  e.g. solve equations, update model, for film, Lagrangian etc.
        virtual void correct();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const IBMForceModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
