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
    nonLinearPlasticModel

Description
    Base class for plastic models with nonLinear-elastic strain

SourceFiles
    nonLinearPlasticModel.C
    nonLinearPlasticModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearPlasticModel_H
#define nonLinearPlasticModel_H

#include "plasticModel.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class nonLinearPlastic Declaration
\*---------------------------------------------------------------------------*/

class nonLinearPlasticModel
:
    public plasticModel
{
    // Private data

        //- Elastic left Cauchy-Green trial strain tensor with volumetric term
        //  removed. Trial means that plasticity is neglected.
        defineVolSurfaceFieldTypeFuncs
        (
            bEbarTrial,
            symmTensor,
            dimensionedSymmTensor("bEbarTrial", dimless, Zero)
        );

        //- Elastic left Cauchy-Green strain tensor with volumetric term removed
        defineVolSurfaceFieldTypeFuncs
        (
            bEbar,
            symmTensor,
            dimensionedSymmTensor("bEbar", dimless, Zero)
        );

        //- Update bEbar consistently with the assumption that det(bEbar) == 1
        //  defaults to off
        const Switch updateBEbarConsistent_;


    // Private Member Functions

        //- Calcualte Ibar such that det(bEbar) == 1
        template<template<class> class PatchField, class Mesh>
        tmp<GeometricField<scalar, PatchField, Mesh>> Ibar
        (
            const GeometricField<symmTensor, PatchField, Mesh>& devBEbar
        );


public:

    // Static data members


    // Constructors

        //- Construct from dictionary
        nonLinearPlasticModel
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~nonLinearPlasticModel();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return the implicit stiffness for a patch
        //  This is the diffusivity for the Laplacian term
        virtual tmp<scalarField> impK(const label) const;

        //- Update the stress
        virtual void correct(volSymmTensorField& sigma)
        {
            correct(sigma, true);
        }

        virtual void correct
        (
            volSymmTensorField& sigma,
            const bool needUpdate
        );

        //- Update the stress surface field
        virtual void correct(surfaceSymmTensorField& sigma)
        {
            correct(sigma, true);
        }

        virtual void correct
        (
            surfaceSymmTensorField& sigma,
            const bool needUpdate
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nonLinearPlasticModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
