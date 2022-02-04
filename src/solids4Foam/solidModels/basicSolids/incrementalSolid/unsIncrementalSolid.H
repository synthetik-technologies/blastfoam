/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    unsIncrementalSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, where the increment of displacement is the
    primary unknown.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    unsIncrementalSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsIncrementalSolid_H
#define unsIncrementalSolid_H

#include "incrementalSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linGeomSolid Declaration
\*---------------------------------------------------------------------------*/

class unsIncrementalSolid
:
    public incrementalSolid
{
protected:
    // Protected data

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Surface displacement gradient
        surfaceTensorField gradDf_;

        //- Surface displacement increment gradient
        surfaceTensorField gradDDf_;


    // Protected member functions

        //- Update displacement fields
        virtual void updateDisplacement();


public:

    // Constructors

        //- Construct from components
        unsIncrementalSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonLinear,
            const bool isSolid = true
        );


    // Destructor

        virtual ~unsIncrementalSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether gradD or gradDD
            //  is the primary solution variable
            virtual surfaceTensorField& solutionGradDf()
            {
                // This model solves for DD
                return gradDDf_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
