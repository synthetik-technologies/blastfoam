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
    Foam::materialModels::elasticPlastic

Description
    Base class for a elasticPlastic material

SourceFiles
    elasticPlasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef elasticPlasticMaterialModel_H
#define elasticPlasticMaterialModel_H

#include "elasticMaterialModel.H"
#include "plasticSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class elasticPlastic Declaration
\*---------------------------------------------------------------------------*/

class elasticPlastic
:
    public elastic,
    public plasticSolver
{
protected:
// Protected data

    // Linear fields

        //- Yield stress
        pointScalarField sigmaY_;
        pointScalarField DsigmaY_;

        //- Equivilant elasticPlastic strain
        pointScalarField epsilonPEq_;
        pointScalarField DepsilonPEq_;

        //- Plastic strain (symm)
        pointSymmTensorField epsilonP_;
        pointSymmTensorField DepsilonP_;


public:

    //- Runtime type information
    TypeName("elasticPlastic")

    // Constructors

        //- Construct from order, size and element type
        elasticPlastic
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const pointVectorField& D,
            const pointVectorField& U,
            const bool planeStress,
            const GeoType geoType
        );

    //- Destructor
    virtual ~elasticPlastic();


    // Member functions

        // Evaluation

            //- Pre update step
            virtual void preUpdate();

            //- Post update step
            virtual void postUpdate(const scalar f);

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
