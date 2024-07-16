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
    Foam::materialModels::elastic

Description
    Base class for a elastic material

SourceFiles
    elasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef elasticMaterialModel_H
#define elasticMaterialModel_H

#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class elastic Declaration
\*---------------------------------------------------------------------------*/

class elastic
:
    public materialModel
{
protected:
// Protected data

    //- Young's modulus
    dimensionedScalar E_;

    //- Poisson's ratio
    dimensionedScalar nu_;

    //- Bulk modulus
    dimensionedScalar K_;

    //- Shear modulus
    dimensionedScalar mu_;

    //- Lame's first parameter
    dimensionedScalar lambda_;


public:

    //- Runtime type information
    TypeName("elastic")

    // Constructors

        //- Construct from order, size and element type
        elastic
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const pointVectorField& D,
            const pointVectorField& U,
            const bool planeStress,
            const GeoType geoType
        );

    //- Destructor
    virtual ~elastic();


    // Member functions

        // Access

            //- Return Young's modulus
            inline const dimensionedScalar& E() const
            {
                return E_;
            }

            //- Return Poisson's ratio
            inline const dimensionedScalar& nu() const
            {
                return nu_;
            }

            //- Return bulk modulu
            inline const dimensionedScalar& K() const
            {
                return K_;
            }

            //- Return shear modulus
            inline const dimensionedScalar& mu() const
            {
                return mu_;
            }

            //- Return Lame's first parameter
            inline const dimensionedScalar& lambda() const
            {
                return lambda_;
            }


        // Evaluation

            //- Compute cauchy stress and return wave speed
            virtual scalar calcSigma
            (
                symmTensor& sigma,
                const label elemi,
                const label rulei
            );

            //- Compute Piola stress and return wave speed
            virtual scalar calcPiola
            (
                tensor& Piola,
                symmTensor& sigma,
                const label elemi,
                const label rulei
            );

            //- Return base wavespeed
            virtual scalar vp() const
            {
                return sqrt
                (
                    (lambda_.value() + 2.0*mu_.value())/rho_.value()
                );
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
