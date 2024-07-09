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
    Foam::materialModels::elasticShell

Description
    Base class for a elastic shell material

SourceFiles
    elasticShellMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef elasticShellMaterialModel_H
#define elasticShellMaterialModel_H

#include "shellMaterialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class elastic Declaration
\*---------------------------------------------------------------------------*/

class elasticShell
:
    public shellMaterialModel
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
    TypeName("elasticShell")

    // Constructors

        //- Construct from order, size and element type
        elasticShell
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const Field<vector>& D,
            const Field<vector>& U,
            const Field<vector>& theta,
            const Field<vector>& omega,
            const GeoType geoType
        );

    //- Destructor
    virtual ~elasticShell();


    // Member functions


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
