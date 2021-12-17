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
    explicitTotalLagSolid

Description:
    Modified version of the explicit solid dynamics library.

    Modifications were made to allow for the use of more complicated
    boundary conditions enforced through displacement, and additional
    generalization


Reference:

    Haider, J., Lee, C. H., Gill, A. J., Bonet, J. (2017),
    "A first-order hyperbolic framework for large strain computational
    solid dynamics: An upwind cell centred Total Lagrangian scheme",
    International Journal for Numerical Methods in Engineering,
    109, 407-456.

    https://github.com/jibranhaider/explicitSolidDynamics

SourceFiles
    explicitTotalLagSolid.C

Also see:
    explicitLinGeomTotalDispSolid.C
    explicitNonLinGeomTotalDispSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitTotalLagSolid_H
#define explicitTotalLagSolid_H

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
                        Class explicitTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitTotalLagSolid
:
    public solidModel
{
    // Private data

        //- Total deformation gradient
        volTensorField F_;

        //- Operations
        operations ops_;

        //- Stabilisation parameter
        scalar beta_;

        //- Is angularMomentum conserved
        Switch angularMomentumConservation_;

        //- Cell displacment
        volVectorField& D_;

        //- Material coordinates
        volVectorField x_;

        //- Face displacement
        surfaceVectorField Df_;

        //- Point displacement
        pointVectorField& pointD_;

        //- Cell velocities
        volVectorField& U_;

        //- linear momentum
        volVectorField rhoU_;

        //- Contact face linear momentum
        surfaceVectorField rhoUC_;

        //- Nodal linear momentum
        pointVectorField pointRhoU_;

        //- Contact traction values
        surfaceVectorField tractionC_;

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

        // Scale factor for Jameson-Schmidt-Turkel smoothing term
        const scalar JSTScaleFactor_;

        //- Energy calculation
        mechanicalEnergies energies_;

        //- Stored impK
        volScalarField impK_;



    // Private Member Functions

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
        explicitTotalLagSolid
        (
            const explicitTotalLagSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const explicitTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }

        //- Does this model take an incremental approach?
        //  i.e. does it solve for DD instead of D?
        //  This defaults to false but can be overwritten
        virtual bool incremental() const
        {
            return false;
        }

public:

    //- Runtime type information
    TypeName("explicitNonLinearGeometryTotalLagrangian");

    // Constructors

        //- Construct from components
        explicitTotalLagSolid
        (
            const word& type,
            dynamicFvMesh& mesh,
            const nonLinearGeometry::nonLinearType&
        );

    // Destructor

        virtual ~explicitTotalLagSolid()
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
