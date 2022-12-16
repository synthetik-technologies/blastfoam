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
    unsLinearTotalDisplacementSolid

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed.

    The stress is calculated by the run-time selectable mechanical law.

    "uns" in the class name refers to the "unstructured" and indicates that
    the face tangential gradient are calculated using a face-Gauss approach.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    unsLinearTotalDisplacementSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsLinearTotalDisplacementSolid_H
#define unsLinearTotalDisplacementSolid_H

#include "UnsLinearGeomSolid.H"
#include "unsTotalDisplacementSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                Class unsLinearTotalDisplacementSolid Declaration
\*---------------------------------------------------------------------------*/

class unsLinearTotalDisplacementSolid
:
    public UnsLinearGeomSolid<unsTotalDisplacementSolid>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        unsLinearTotalDisplacementSolid(const unsLinearTotalDisplacementSolid&);

        //- Disallow default bitwise assignment
        void operator=(const unsLinearTotalDisplacementSolid&);


public:

    //- Runtime type information
    TypeName("unsLinearTotalDisplacement");

    // Constructors

        //- Construct from components
        unsLinearTotalDisplacementSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~unsLinearTotalDisplacementSolid()
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
