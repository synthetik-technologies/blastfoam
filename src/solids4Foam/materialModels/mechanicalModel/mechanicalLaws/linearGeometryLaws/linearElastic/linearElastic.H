/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

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
    linearElastic

Description
    Hookean linear elastic mechanical law.

SourceFiles
    linearElastic.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef linearElastic_H
#define linearElastic_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class linearElastic
:
    public mechanicalLaw
{
    // Private data

        //- Shear modulus (second Lame parameter)
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- First Lame parameter
        dimensionedScalar lambda_;

        //- Total strain (vol)
        defineVolSurfaceFieldFuncs
        (
            epsilon,
            symmTensor,
            dimensionedSymmTensor("zero", dimless, symmTensor::zero)
        );

        //- Initial (residual) stress field
        //  Defaults to zero if not specified
        defineVolSurfaceFieldFuncs
        (
            sigma0,
            symmTensor,
            dimensionedSymmTensor("zero", dimPressure, symmTensor::zero)
        );

    // Private Member Functions

        //- Calculate hydrostatic component of the stress tensor
        void calculateHydrostaticStress
        (
            volScalarField& sigmaHyd,
            const volScalarField& trEpsilon
        );

        //- Calculate hydrostatic component of the stress tensor
        void calculateHydrostaticStress
        (
            surfaceScalarField& sigmaHyd,
            const surfaceScalarField& trEpsilon
        );

        //- Disallow default bitwise copy construct
        linearElastic(const linearElastic&);

        //- Disallow default bitwise assignment
        void operator=(const linearElastic&);

public:

    //- Runtime type information
    TypeName("linearElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        linearElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~linearElastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the implicit stiffness for a patch
        //  This is the diffusivity for the Laplacian term
        virtual tmp<scalarField> impK(const label) const;

        //- Return bulk modulus
        virtual tmp<volScalarField> bulkModulus() const;

        //- Return elastic modulus
        virtual tmp<volScalarField> elasticModulus() const;

        //- Return shear modulus
        virtual tmp<volScalarField> shearModulus() const;

        //- Provide access to mu for the coupledUnsNonLinLinearElasticSolid
        virtual const dimensionedScalar& mu() const;

        //- Bulk modulus
        virtual const dimensionedScalar& K() const;

        //- Young's modulus
        virtual const dimensionedScalar& E() const;

        //- Poisson's ratio
        virtual const dimensionedScalar& nu() const;

        //- Lambda
        virtual const dimensionedScalar& lambda() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
