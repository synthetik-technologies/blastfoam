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
    explicitRiemannSolid

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
    explicitRiemannSolid.C

Also see:
    explicitLinRiemannSolid.C
    explicitNonLinRiemannSolid.C

\*---------------------------------------------------------------------------*/

#ifndef explicitRiemannSolid_H
#define explicitRiemannSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "mechanicalEnergies.H"
#include "operations.H"
#include "mechanics.H"
#include "gradientSchemes.H"
#include "angularMomentum.H"
#include "dynamicRelaxation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class explicitRiemannSolid Declaration
\*---------------------------------------------------------------------------*/

class explicitRiemannSolid
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

        //- Increment of cell displacment
        volVectorField& DD_;

        //- Gradient of cell displacment
        volTensorField& gradD_;

        //- Gradient of the increment of cell displacment
        volTensorField& gradDD_;

        //- Material coordinates
        volVectorField x_;
        surfaceVectorField xf_;
        pointVectorField xN_;

        //- Point displacement
        pointVectorField& pointD_;

        //- Increment of point displacement
        pointVectorField& pointDD_;

        //- Density
        const volScalarField& rho_;

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

        //- Angular momentum
        angularMomentum am_;

        //- Dynamic relaxation
        dynamicRelaxation relaxation_;

        //- Piola stress tensor
        volTensorField P_;

        //- Compressive wave speed field
        volScalarField pWaveSpeed_;

        //- Shear wave speed field
        volScalarField sWaveSpeed_;

        //- Include stabilisation terms
        bool useStabilisation_;

        //- Include bulk viscosity
        bool useBulkViscosity_;

        //- Energy calculation
        mechanicalEnergies energies_;

        //- Use the penalised-TOUCH correction
        bool pTouch_;

        //- Use filtering (aka solve nodal linear momentum)
        bool filter_;

        //- Stored impK
        volScalarField impK_;
        surfaceScalarField impKf_;

        //- Current time index
        label curIndex_;



    // Private Member Functions

        //- Smooth the hydrostatic pressure field
        void solveGEqns
        (
            volVectorField& rhoURHS,
            const label stage
        );

        //- Return smoothed U old time based on Lax-Friedlichs method
        void updateFluxes();

        //- Disallow default bitwise copy construct
        explicitRiemannSolid
        (
            const explicitRiemannSolid&
        );

        //- Disallow default bitwise assignment
        void operator=(const explicitRiemannSolid&);


public:

    //- Runtime type information
    ClassName("explicitRiemann");

    // Constructors

        //- Construct from components
        explicitRiemannSolid
        (
            const word& type,
            fvMesh& mesh,
            const nonLinearGeometry::nonLinearType&
        );

    // Destructor

        virtual ~explicitRiemannSolid()
        {}


    // Member Functions

        // Access

            //- Does this model take an incremental approach?
            //  i.e. does it solve for DD instead of D?
            //  This defaults to false but can be overwritten
            virtual bool incremental() const
            {
                return true;
            }

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual const volVectorField& solutionD() const
            {
                return DD();
            }

            //- Each solidModel must indicate whether gradD or gradDD
            //  is the primary solution variable
            virtual const volTensorField& solutionGradD() const
            {
                return gradDD();
            }

            //- Return the patch implicit stiffness
            virtual tmp<scalarField> impK(const fvPatch& patch) const
            {
                return impKf_.boundaryField()[patch.index()];
            }

        // Edit

            //- Decode
            void decode();

            //- Update fields
            void update(const bool correctSigma = true);

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            virtual tmp<Foam::vectorField> nf(const fvPatch& patch) const
            {
                return mech_.n().boundaryField()[patch.index()];
            }

            //- Return the current Courant number
            virtual scalar CoNum() const;

            //- Return the max Courant number
            virtual scalar maxCoNum() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
