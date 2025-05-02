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
    Foam::materialModel

Description
    Base class for a finite element

SourceFiles
    materialModelI.H
    materialModel.C

\*---------------------------------------------------------------------------*/

#ifndef materialModel_H
#define materialModel_H

#include "feMesh1.H"
#include "scalarMatrices.H"
#include "syncTools.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class materialModel Declaration
\*---------------------------------------------------------------------------*/

class materialModel
{
public:

    enum GeoType
    {
        LINEAR,
        UPDATE_LAGRANGIAN,
        TOTAL_LAGRANGIAN
    };

protected:

    //- Reference to the feMesh1
    const feMesh1& mesh_;

    //- Geometry type
    GeoType geoType_;

    //- Displacement field
    const pointVectorField& D_;

    //- Velocity field
    const pointVectorField& U_;

    //- Density
    dimensionedScalar rho_;

    //- Name of the material
    word name_;

    //- Plane stress
    bool planeStress_;

    //- Plane stress direction
    direction planeStressDir_;

    //- Non-plane stress directions
    direction nonPlaneStressDir1_;
    direction nonPlaneStressDir2_;

    // Update variables
    mutable bool update_;


public:

    //- Runtime type information
    TypeName("materialModel")

    // Declare runtime constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            materialModel,
            linear,
            (
                const dictionary& dict,
                const feMesh1& mesh,
                const pointVectorField& D,
                const pointVectorField& U,
                const bool planeStress,
                const GeoType geoType
            ),
            (dict, mesh, D, U, planeStress, geoType)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            materialModel,
            nonLinear,
            (
                const dictionary& dict,
                const feMesh1& mesh,
                const pointVectorField& D,
                const pointVectorField& U,
                const bool planeStress,
                const GeoType geoType
            ),
            (dict, mesh, D, U, planeStress, geoType)
        );

    // Constructors

        //- Construct from order, size and element type
        materialModel
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const pointVectorField& D,
            const pointVectorField& U,
            const bool planeStress,
            const GeoType geoType
        );


    // Selectors



        //- Select Update-Lagrangian materialModel
        static autoPtr<materialModel> New
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const pointVectorField& D,
            const pointVectorField& U,
            const bool planeStress,
            const GeoType
        );

//         //- Select Update-Lagrangian materialModel
//         static autoPtr<materialModel> New
//         (
//             const dictionary& dict,
//             const feMesh1& mesh,
//             const pointVectorField& D,
//             const pointVectorField& U,
//             const bool planeStress
//         );


    //- Destructor
    virtual ~materialModel();


    // Member functions

        // Access

            //- Return the geometry type
            inline GeoType geoType() const
            {
                return geoType_;
            }

            //- Return the material name
            inline const word& name() const
            {
                return name_;
            }

            //- Return the density
            inline const dimensionedScalar& rho() const
            {
                return rho_;
            }

            //- Return elements that use this material
            virtual const labelList& elements() const;


        // Evaluation

            //- Pre update step
            virtual void preUpdate() {}

            //- Post update step
            virtual void postUpdate(const scalar f) {}

            //- Compute the reference gradient
            template<template<class> class ListT>
            static vector refGradient
            (
                const ListT<scalar>& f,
                const scalarRectangularMatrix& dshape
            )
            {
                // Compute gradient in reference space
                vector gradf(Zero);
                for (label dimi = 0; dimi < dshape.n(); dimi++)
                {
                    scalar gf = 0.0;
                    for (label si = 0; si < dshape.m(); si++)
                    {
                        gf += f[si]*dshape(si, dimi);
                    }
                    gradf[dimi] = gf;
                }
                return gradf;
            }

            //- Compute the physical gradient
            template<template<class> class ListT>
            static vector gradient
            (
                const ListT<scalar>& f,
                const scalarRectangularMatrix& dshape,
                const tensor& invJ
            )
            {
                return refGradient(f, dshape) & invJ;
            }

            //- Compute the reference gradient
            template<template<class> class ListT>
            static tensor refGradient
            (
                const ListT<vector>& f,
                const scalarRectangularMatrix& dshape
            )
            {
                // Compute gradient in reference space
                tensor gradf(Zero);
                for (label dimi = 0; dimi < 3; dimi++)
                {
                    for (label dimj = 0; dimj < dshape.n(); dimj++)
                    {
                        scalar gf = 0.0;
                        for (label si = 0; si < dshape.m(); si++)
                        {
                            gf += f[si][dimi]*dshape(si, dimj);
                        }
                        gradf(dimi, dimj) = gf;
                    }
                }
                return gradf;
            }

            //- Compute the physical gradient
            template<template<class> class ListT>
            static tensor gradient
            (
                const ListT<vector>& f,
                const scalarRectangularMatrix& dshape,
                const tensor& invJ
            )
            {
                return refGradient(f, dshape) & invJ;
            }

            //- Compute strain tensor
            void updateStrain
            (
                symmTensor& epsilon,
                const scalar coeff,
                const symmTensor& sigma
            ) const;

            //- Compute the deformation gradient tensor
            template<template<class> class ListT>
            static tensor calcF
            (
                const ListT<vector>& D,
                const scalarRectangularMatrix& dshape,
                const tensor& invJ
            )
            {
                // Transform gradient to physical space
                return gradient(D, dshape, invJ) + tensor::I;
            }

            //- Compute cauchy stress and return wave speed
            virtual scalar calcSigma
            (
                symmTensor& sigma,
                const label elemi,
                const label rulei
            ) = 0;

            //- Compute Piola stress and return wave speed
            virtual scalar calcPiola
            (
                tensor& Piola,
                symmTensor& sigma,
                const label elemi,
                const label rulei
            ) = 0;

            //- Add force, return wave speed
            virtual scalar addForce
            (
                UList<vector>& f,
                UList<symmTensor>& sigma,
                const label elemi
            );

            //- Return base wavespeed
            virtual scalar vp() const = 0;



    // Static member functions

        //- Compute bulk and shear moduli from Young's modulus and
        //  Poisson ratio
        static void Kmu_Enu
        (
            const scalar E, const scalar nu,
            scalar& K, scalar& mu,
            const bool ps
        );

        //- Compute Young's modulus and Poisson ratio from bulk and
        //  shear moduli
        static void Enu_Kmu
        (
            const scalar K, const scalar mu,
            scalar& E, scalar& nu,
            const bool ps
        );

        //- Compute Lame's first parameter from bulk and shear modulus
        static scalar lambda_Kmu
        (
            const scalar K, const scalar mu,
            const bool ps
        );

        //- Compute Lame's first parameter from Young's modulus and
        //  Poisson ratio
        static scalar lambda_Enu
        (
            const scalar E, const scalar nu,
            const bool ps
        );

        virtual void write() const
        {}

        template<class Type>
        void writeIpField
        (
            const word& name,
            const dimensionSet& dims,
            const List<Type> ipField
        ) const
        {
            GeometricField<Type, pointPatchField, pointMesh> pField
            (
                IOobject
                (
                    name,
                    mesh_.time().name(),
                    mesh_.mesh()
                ),
                mesh_.pMesh(),
                dimensioned<Type>(name, dims, Zero)
            );


            forAll(mesh_.elements(), elemi)
            {
                const element& elem = mesh_.elements()[elemi];
                const labelList& ips = mesh_.ipLabels()[elemi];
                const List<scalarList>& shapes = mesh_.shapes()[elemi];
                const List<scalar>& Ws = mesh_.Ws()[elemi];
                UIndirectList<Type> pField_loc(pField, elem);

                forAll(elem.ir(), rulei)
                {
                    const scalar w = elem.ir()[rulei].w()*Ws[rulei];
                    const scalarList& shape = shapes[rulei];
                    forAll(shape, si)
                    {
                        pField_loc[si] += w*ipField[ips[rulei]]*shape[si];
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh_.mesh(),
                pField,
                plusEqOp<Type>(),
                pTraits<Type>::zero
            );

            pField.primitiveFieldRef() /= mesh_.W();
            mesh_.pushUntransformedData(pField);

            pField.write();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
