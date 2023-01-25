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
    linearTotalDisplacementSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement is the primary
    unknown.

    The stress is calculated by the run-time selectable mechanical law.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    linearTotalDisplacementSolid.C

\*---------------------------------------------------------------------------*/

#ifndef linearTotalDisplacementSolid_H
#define linearTotalDisplacementSolid_H

#include "LinearGeomSolid.H"
#include "totalDisplacementSolids.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class linearTotalDisplacementSolid Declaration
\*---------------------------------------------------------------------------*/

class linearTotalDisplacementSolid
:
    public LinearGeomSolid<totalDisplacementSolid>
{
    // Private data

        //- Predict new time-step fields using the velocity field
        const Switch predictor_;


    // Private Member Functions

        //- Predict the fields for the next time-step based on the
        //  previous time-steps
        void predict();

        //- Disallow default bitwise copy construct
        linearTotalDisplacementSolid(const linearTotalDisplacementSolid&);

        //- Disallow default bitwise assignment
        void operator=(const linearTotalDisplacementSolid&);


public:

    //- Runtime type information
    TypeName("linearTotalDisplacement");

    // Constructors

        //- Construct from components
        linearTotalDisplacementSolid(dynamicFvMesh& mesh);

        //- Construct from components and typename
        linearTotalDisplacementSolid(const word& type, dynamicFvMesh& mesh);

    // Destructor

        virtual ~linearTotalDisplacementSolid()
        {}


    // Member Functions

        //- Evolve the solid solver and solve the mathematical model
        virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
