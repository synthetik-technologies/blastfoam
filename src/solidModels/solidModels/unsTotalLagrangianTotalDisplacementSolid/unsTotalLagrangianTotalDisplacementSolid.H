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
    unsTotalLagrangianTotalDisplacementSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and total displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    "uns" in the class name refers to the "unstructured" and indicates that
    the face tangential gradient are calculated using a face-Gauss approach.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    unsTotalLagrangianTotalDisplacementSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsTotalLagrangianTotalDisplacementSolid_H
#define unsTotalLagrangianTotalDisplacementSolid_H

#include "UnsTotalLagrangianGeomSolid.H"
#include "totalDisplacementSolids.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
        Class unsTotalLagrangianTotalDisplacementSolid Declaration
\*---------------------------------------------------------------------------*/

class unsTotalLagrangianTotalDisplacementSolid
:
    public UnsTotalLagrangianGeomSolid<unsTotalDisplacementSolid>
{
    // Private data

        //- Switch to enable/disable nonlinear geometry
        const Switch nonLinear_;

        //- Viscosity/damping coefficient
        const dimensionedScalar K_;

        //- Relative tolerance
        scalar relativeTol_;

    // Private Member Functions

        //- Return the displacement residual
        scalar residual(const volVectorField& D) const;

        //- Disallow default bitwise copy construct
        unsTotalLagrangianTotalDisplacementSolid
        (
            const unsTotalLagrangianTotalDisplacementSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const unsTotalLagrangianTotalDisplacementSolid&);

        virtual void readDict();


public:

    //- Runtime type information
    TypeName("unsTotalLagrangianTotalDisplacement");

    // Constructors

        //- Construct from components
        unsTotalLagrangianTotalDisplacementSolid(fvMesh& mesh);


    // Destructor

        virtual ~unsTotalLagrangianTotalDisplacementSolid()
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
