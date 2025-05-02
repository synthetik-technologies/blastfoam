/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::shellMaterialModel

Description
    Base class for a shell material model

SourceFiles
    shellMaterialModelI.H
    shellMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef shellMaterialModel_H
#define shellMaterialModel_H

#include "feMesh1.H"
#include "materialModel.H"
#include "scalarMatrices.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class shellMaterialModel Declaration
\*---------------------------------------------------------------------------*/

class shellMaterialModel
:
    public materialModel
{
protected:

    //- Rotation angle field
    const Field<vector>& theta_;

    //- Angular velocity field
    const Field<vector>& omega_;

    //- Thickness
    const scalar h_;


public:

    //- Runtime type information
    TypeName("shellMaterialModel")

    // Declare runtime constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            shellMaterialModel,
            linear,
            (
                const dictionary& dict,
                const feMesh1& mesh,
                const Field<vector>& D,
                const Field<vector>& U,
                const Field<vector>& theta,
                const Field<vector>& omega,
                const GeoType geoType
            ),
            (dict, mesh, D, U, theta, omega, geoType)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            shellMaterialModel,
            nonLinear,
            (
                const dictionary& dict,
                const feMesh1& mesh,
                const Field<vector>& D,
                const Field<vector>& U,
                const Field<vector>& theta,
                const Field<vector>& omega,
                const GeoType geoType
            ),
            (dict, mesh, D, U, theta, omega, geoType)
        );

    // Constructors

        //- Construct from order, size and element type
        shellMaterialModel
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const Field<vector>& D,
            const Field<vector>& U,
            const Field<vector>& theta,
            const Field<vector>& omega,
            const GeoType geoType
        );


    // Selectors



        //- Select Update-Lagrangian shellMaterialModel
        static autoPtr<shellMaterialModel> New
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const Field<vector>& D,
            const Field<vector>& U,
            const Field<vector>& theta,
            const Field<vector>& omega,
            const GeoType
        );

    //- Destructor
    virtual ~shellMaterialModel();


    // Member functions

        //- Return elements that use this material
        virtual const labelList& elements() const;

        // Evaluation

            //- Pre update step
            virtual void preUpdate() {}

            //- Post update step
            virtual void postUpdate() {}

            //- Compute strain tensor
            void updateStrain
            (
                symmTensor& epsilon,
                const scalar coeff,
                const symmTensor& sigma
            ) const;

            //- Compute covariant basis vectors
            FixedList<vector, 3> covariantBasis
            (
                const tensor& gradX,
                const vector& v
            ) const;

            //- Compute contravariant basis vectors
            FixedList<vector, 3> contravariantBasis
            (
                const FixedList<vector, 3>& cov
            ) const;

            //- Compute local basis vectors
            FixedList<vector, 3> localBasis
            (
                const FixedList<vector, 3>& cov
            ) const;

            //- Compute the deformation gradient tensor
            template<template<class> class ListT>
            tensor calcF
            (
                const ListT<vector>& D,
                const ListT<vector>& theta,
                const vector& n,
                const scalar h,
                const scalarList& shape,
                const scalarRectangularMatrix& dshape,
                const tensor& invJ
            ) const
            {
                vector th(Zero);
                forAll(shape, i)
                {
                    th += theta[i]*shape[i];
                }

                // Transform gradient to physical space
                return
                    gradient(D, dshape, invJ)
                  + h*gradient(theta, dshape, invJ)
                  + th*n
                  + tensor::I;
            }

            //- Add force and moment, return wave speed
            virtual scalar addForce
            (
                UList<vector>& f,
                UList<vector>& moment,
                UList<symmTensor>& sigma,
                const label elemi
            );

            virtual scalar addForce
            (
                UList<vector>& f,
                UList<symmTensor>& sigma,
                const label elemi
            )
            {
                NotImplemented;
                return 0;
            }

            //- Return thickness
            scalar thickness() const
            {
                return h_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
