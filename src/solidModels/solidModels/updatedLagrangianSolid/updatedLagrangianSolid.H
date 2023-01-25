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
    updatedLagrangianSolid

Description
    Mathematical model for finite strains and rotations, where an updated
    Lagrangian and incremental displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    More details in:
    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    updatedLagrangianSolid.C

\*---------------------------------------------------------------------------*/

#ifndef updatedLagrangianSolid_H
#define updatedLagrangianSolid_H

#include "UpdatedLagrangianGeomSolid.H"
#include "incrementalSolids.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                    Class updatedLagrangianSolid Declaration
\*---------------------------------------------------------------------------*/

class updatedLagrangianSolid
:
    public UpdatedLagrangianGeomSolid<incrementalSolid>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        updatedLagrangianSolid(const updatedLagrangianSolid&);

        //- Disallow default bitwise assignment
        void operator=(const updatedLagrangianSolid&);


public:

    //- Runtime type information
    TypeName("updatedLagrangian");

    // Constructors

        //- Construct from components
        updatedLagrangianSolid(dynamicFvMesh& mesh);


    // Destructor

        virtual ~updatedLagrangianSolid()
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
