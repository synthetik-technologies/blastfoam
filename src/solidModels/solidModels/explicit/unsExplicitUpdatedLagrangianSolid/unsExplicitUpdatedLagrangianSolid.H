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
    unsExplicitUpdatedLagrangianSolid

Description
    Version of the unsLinGeomSolid solid model, where a fully explicit
    solution algorithm is used.

    The explicit implementation uses central differencing, where a linear bulk
    viscosity pressure term dissipates high frequencies.

    Details of the central differencing solution algorithm can be found in, e.g.
    - Belytschko, Liu, Moran, 2000, Box 6.1
    - Ding and Ye, Laser shock peening, Performance and process simulation,
      2006, Section 3.5.2

    The linear bulk viscosity term is described in the Abaqus documentation.

    A Jameson-Schmidt-Turkel (JST) 4th order diffusion term is used for
    stabilisation.

Author
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    unsExplicitUpdatedLagrangianSolid.C

\*---------------------------------------------------------------------------*/

#ifndef unsExplicitUpdatedLagrangianSolid_H
#define unsExplicitUpdatedLagrangianSolid_H

#include "UnsUpdatedLagrangianGeomSolid.H"
#include "unsExplicitNonLinearSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{

/*---------------------------------------------------------------------------*\
                Class unsExplicitUpdatedLagrangianSolid Declaration
\*---------------------------------------------------------------------------*/

class unsExplicitUpdatedLagrangianSolid
:
    public UnsUpdatedLagrangianGeomSolid<unsExplicitNonLinearSolid>
{
    // Private Member Functions

        //- Return the face traction
        virtual tmp<surfaceVectorField> tractionSf() const;

        //- Return the divergence of stress
        virtual tmp<volVectorField> divStress() const;

        //- Disallow default bitwise copy construct
        unsExplicitUpdatedLagrangianSolid
        (
            const unsExplicitUpdatedLagrangianSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const unsExplicitUpdatedLagrangianSolid&);


public:

    //- Runtime type information
    TypeName("unsExplicitUpdatedLagrangian");

    // Constructors

        //- Construct from components
        unsExplicitUpdatedLagrangianSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~unsExplicitUpdatedLagrangianSolid()
        {}


    // Member Functions

        // Edit

            //- Update stress
            virtual void update(const bool correctSigma = true);

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
