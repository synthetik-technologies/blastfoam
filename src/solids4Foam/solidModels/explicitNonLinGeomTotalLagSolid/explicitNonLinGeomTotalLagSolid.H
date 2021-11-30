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
    explicitNonLinGeomTotalLagSolid

Description
    Version of the linGeomTotalDispSolid solid model, where a fully explicit
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
    explicitNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitNonLinGeomTotalLagSolid_H
#define explicitNonLinGeomTotalLagSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "mechanicalEnergies.H"
#include "operations.H"
#include "mechanics.H"
#include "gradientSchemes.H"
#include "interpolationSchemes.H"
#include "angularMomentum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class explicitNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitNonLinGeomTotalLagSolid
:
    public solidModel
{
    // Private data

        //- Total deformation gradient
        volTensorField F_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Operations
        operations ops_;

        //- Stabilisation parameter
        scalar beta_;

        //- Is angularMomentum conserved
        Switch angularMomentumConservation_;

        //- Spatial cell centres
        volVectorField x_;

        //- Spatial face centres
        surfaceVectorField xf_;

        //- Material nodal locations
        pointVectorField pointX_;

        //- Face velocities
        surfaceVectorField Uf_;

        //- linear momentum
        volVectorField rhoU_;

        //- Face linear momentum
        surfaceVectorField rhoUf_;

        //- Nodal linear momentum
        pointVectorField pointRhoU_;

        //- Cell centred traction
        volVectorField traction_;

        //- Wave speed field
        surfaceVectorField tractionf_;

        //- Mechanics
        mechanics mech_;

        //- Interpolation schemes
        interpolationSchemes interpSchemes_;

        //- Gradients
        gradientSchemes gradSchemes_;

        //- Angular momentum
        angularMomentum am_;

        //- Piola stress tensor
        volTensorField P_;

        //- Compressive wave speed field
        volScalarField pWaveSpeed_;

        //- Shear wave speed field
        volScalarField sWaveSpeed_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;

        //- Energy calculation
//         mechanicalEnergies energies_;




    // Private Member Functions

        // Linear momentum boundary types
        wordList linearMomentumBoundaryTypes() const;

        // Traction boundary types
        wordList tractionBoundaryTypes() const;

        //- Update the stress field
        void updateStress();

        //- Smooth the hydrostatic pressure field
        void solveGEqns
        (
            volVectorField& rhoURHS,
            volVectorField& rhoURHS1,
            const label stage
        );

        //- Return smoothed U old time based on Lax-Friedlichs method
        void updateFluxes();

        //- Disallow default bitwise copy construct
        explicitNonLinGeomTotalLagSolid
        (
            const explicitNonLinGeomTotalLagSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const explicitNonLinGeomTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("explicitNonLinearGeometryTotalLagrangian");

    // Constructors

        //- Construct from components
        explicitNonLinGeomTotalLagSolid(dynamicFvMesh& mesh);

    // Destructor

        virtual ~explicitNonLinGeomTotalLagSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for D
                return D();
            }

        // Edit

            //- Update the size of the time-step
            virtual void setDeltaT(Time& runTime);

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Return the current Courant number
            virtual scalar CoNum() const;

            //- Return the max Courant number
            virtual scalar maxCoNum() const;

            //- End of the program
            virtual void end()
            {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
