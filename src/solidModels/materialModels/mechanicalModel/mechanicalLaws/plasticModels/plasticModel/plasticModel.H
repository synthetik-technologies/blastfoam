/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
31-03-2022 Synthetik Applied Technologies: Added base plasticModel
-------------------------------------------------------------------------------
License
    This file is a derivative work of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    plasticModel

Description
    Base class for a plastic model

SourceFiles
    plasticModel.C


\*---------------------------------------------------------------------------*/

#ifndef plasticModel_H
#define plasticModel_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class plasticModel
:
    public mechanicalLaw
{
protected:

    // Protected Data

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        virtual dimensionedScalar sigmaY0() const = 0;

        //- Total strain
        defineVolSurfaceFieldTypeFuncs
        (
            epsilon,
            symmTensor,
            dimensionedSymmTensor("epsilon", dimless, Zero)
        );

        //- Yield stress
        defineVolSurfaceFieldTypeFuncs
        (
            sigmaY,
            scalar,
            this->sigmaY0()
        );

        //- Total plastic strain
        defineVolSurfaceFieldTypeFuncs
        (
            epsilonP,
            symmTensor,
            dimensionedSymmTensor("epsilonP", dimless, Zero)
        );

        //- Incremental change of plastic strain
        defineVolSurfaceFieldTypeFuncs
        (
            DEpsilonP,
            symmTensor,
            dimensionedSymmTensor("DEpsilonP", dimless, Zero)
        );

        //- Equivalent plastic strain increment
        defineVolSurfaceFieldTypeFuncs
        (
            DEpsilonPEq,
            scalar,
            dimensionedScalar("DEpsilonPEq", dimless, Zero)
        );

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        defineVolSurfaceFieldTypeFuncs
        (
            DLambda,
            scalar,
            dimensionedScalar("DLambda", dimless, Zero)
        );

        //- Equivalent plastic strain
        defineVolSurfaceFieldTypeFuncs
        (
            epsilonPEq,
            scalar,
            dimensionedScalar("epsilonPEq", dimless, Zero)
        );

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- plasticN is the return direction to the yield surface
        defineVolSurfaceFieldTypeFuncs
        (
            plasticN,
            symmTensor,
            dimensionedSymmTensor("plasticN", dimless, Zero)
        );

        //- Linear plastic modulus. It is only used if plasticity is linear,
        //  defined by two points on the stress-plastic strain curve
        scalar Hp_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Protected Member Functions

        //- Does the model have nonlinear plasticity
        virtual bool nonLinearPlasticity() const = 0;

        //- Disallow default bitwise copy construct
        plasticModel(const plasticModel&);

        //- Disallow default bitwise assignment
        void operator=(const plasticModel&);


    //- Protected member functions for rate dependent models

        //- Set local values for fields not passed in functions
        virtual void setCellValues(const label celli) const
        {}
        virtual void setVolPatchFaceValues
        (
            const label patchi,
            const label facei
        ) const
        {}
        virtual void setFaceValues(const label facei) const
        {}
        virtual void setSurfacePatchFaceValues
        (
            const label patchi,
            const label facei
        ) const
        {}

        //- Return the current yield stress
        virtual scalar curYieldStress
        (
            const scalar epsilonPEqOld, // Old equivalent plastic strain
            const scalar curEpsilonPEq, // Current equivalent plastic strain
            const scalar J = 1.0        // Current Jacobian
        ) const = 0;

        //- Evaulate current value of the yield function
        virtual scalar yieldFunction
        (
            const scalar epsilonPEqOld, // Old equivalent plastic strain
            const scalar magSTrial,     // Deviatoric trial stress magnitude
            const scalar DLambda,       // Plastic multiplier
            const scalar muBar,         // Scaled shear modulus
            const scalar J = 1.0        // Current Jacobian
        ) const;

        //- Iteratively calculate plastic multiplier increment (DLambda)
        //  and current yield stress using Newton's method
        virtual void newtonLoop
        (
            scalar& DLambda,            // Plastic multiplier
            scalar& curSigmaY,          // Current yield stress
            const scalar epsilonPEqOld, // Old equivalent plastic strain
            const scalar magSTrial,     // Deviatoric trial stress magnitude
            const scalar muBar,         // Scaled shear modulus
            const scalar maxMagDEpsilon,// Max strain increment magnitude
            const scalar J = 1.0        // Current Jacobian
        ) const;

        //- Update plasticN, DLambda, DSigmaY and sigmaY
        virtual void updatePlasticity
        (
            symmTensor& plasticN,       // Plastic return direction
            scalar& DLambda,            // Plastic multiplier increment
            scalar& sigmaY,             // Yield stress
            const scalar sigmaYOld,     // Yield stress old time
            const scalar fTrial,        // Trial yield function
            const symmTensor& sTrial,   // Trial deviatoric stress
            const scalar epsilonPEqOld, // Old equivalent plastic strain
            const scalar muBar,         // Scaled shear modulus
            const scalar maxMagDEpsilon,// Max strain increment magnitude
            const scalar J = 1.0        // Current Jacobian
        ) const;

public:

    //- Runtime type information
    TypeName("plasticModel");

    // Static data members


    // Constructors

        //- Construct from dictionary
        plasticModel
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~plasticModel();


    // Member Functions

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return elastic modulus
        virtual tmp<volScalarField> elasticModulus() const;

        //- Return shear modulus
        virtual tmp<volScalarField> shearModulus() const;

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();

        //- Return the desired new time-step
        virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
