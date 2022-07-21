/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    unsTotalDispSolid

Description
    Base class for unstructured total displacement solid models

SourceFiles
    unsTotalDispSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsTotalDispSolid_H
#define unsTotalDispSolid_H

#include "totalDispSolid.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class unsTotalDispSolid Declaration
\*---------------------------------------------------------------------------*/

class unsTotalDispSolid
:
    public totalDispSolid
{
protected:
    // Protected data

        //- Cauchy (true) stress surface field
        surfaceSymmTensorField sigmaf_;

        //- Surface displacement gradient
        surfaceTensorField gradDf_;

        //- Surface displacement increment gradient
        surfaceTensorField gradDDf_;


    // Protected Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void updateDisplacement();


public:

    // Constructors

        //- Construct from components
        unsTotalDispSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonLinear,
            const bool isSolid = true
        );

    // Destructor

        virtual ~unsTotalDispSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether gradD or gradDD
            //  is the primary solution variable
            virtual const surfaceTensorField& solutionGradDf() const
            {
                // This model solves for DD
                return gradDf_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
