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
    explicitLinearSolid

Description
    Version of the linearSolid solid model, where a fully explicit
    solution algorithm is used.

    The explicit implementation uses central differencing, where a linear bulk
    viscosity pressure term dissipates high frequencies.

    Details of the central differencing solution algorithm can be found in, e.g.
    - Belytschko, Liu, Moran, 2000, Box 6.1
    - Ding and Ye, Laser shock peening, Performance and process simulation,
      2006, Section 3.5.2

    The bulk viscosity term is described in the Abaqus documentation.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    explicitLinearSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitLinearSolid_H
#define explicitLinearSolid_H

#include "explicitSolid.H"
#include "LinearGeomSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class explicitLinearSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitLinearSolid
:
    public LinearGeomSolid<explicitSolid>
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        explicitLinearSolid(const explicitLinearSolid&);

        //- Disallow default bitwise assignment
        void operator=(const explicitLinearSolid&);

protected:

    // Protected member functions

        //- Correct velocity boundary conditions
        void correctUBoundaryConditions(volVectorField& U);

        //- Return the divergence of stress
        virtual tmp<volVectorField> divStress() const;

        //- Return the face traction
        virtual tmp<surfaceVectorField> tractionSf() const;

public:

    //- Runtime type information
    TypeName("explicitLinear");

    // Constructors

        //- Construct from components
        explicitLinearSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~explicitLinearSolid()
        {}


    // Member Functions

        // Access


        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Update the stress field
            virtual void update(const bool correctSigma = true);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
