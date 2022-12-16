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
    unsUpdatedLagrangianSolid

Description
    Mathematical model for finite strains and rotations, where an updated
    Lagrangian and incremental displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    "uns" in the class name refers to the "unstructured" and indicates that
    the face tangential gradient are calculated using a face-Gauss approach.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    unsUpdatedLagrangianSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsUpdatedLagrangianSolid_H
#define unsUpdatedLagrangianSolid_H

#include "unsIncrementalSolid.H"
#include "UnsUpdatedLagrangianGeomSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class unsUpdatedLagrangianSolid Declaration
\*---------------------------------------------------------------------------*/

class unsUpdatedLagrangianSolid
:
    public UnsUpdatedLagrangianGeomSolid<unsIncrementalSolid>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        unsUpdatedLagrangianSolid(const unsUpdatedLagrangianSolid&);

        //- Disallow default bitwise assignment
        void operator=(const unsUpdatedLagrangianSolid&);

public:

    //- Runtime type information
    TypeName("unsUpdatedLagrangian");

    // Constructors

        //- Construct from components
        unsUpdatedLagrangianSolid(dynamicFvMesh& mesh);


    // Destructor

        virtual ~unsUpdatedLagrangianSolid()
        {}


    // Member Functions

    // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
