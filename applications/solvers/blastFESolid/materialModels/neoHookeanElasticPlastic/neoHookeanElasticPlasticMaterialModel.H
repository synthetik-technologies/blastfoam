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
    Foam::materialModels::neoHookeanElasticPlastic

Description
    Base class for a neoHookeanElasticPlastic material

SourceFiles
    neoHookeanElasticPlasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElasticPlasticMaterialModel_H
#define neoHookeanElasticPlasticMaterialModel_H

#include "elasticPlasticMaterialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class neoHookeanElasticPlastic Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElasticPlastic
:
    public elasticPlastic
{
protected:
// Protected data

    //- bBar
    pointSymmTensorField bBar_;
    pointSymmTensorField DbBar_;

    //- Consistently update bBar
    bool bBarConsistent_;


    // Protected member functions

        //- Calulcate IBar
        scalar calcIBar(const symmTensor& EBar) const;


public:

    //- Runtime type information
    TypeName("neoHookeanElasticPlastic")

    // Constructors

        //- Construct from order, size and element type
        neoHookeanElasticPlastic
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const pointVectorField& D,
            const pointVectorField& U,
            const bool planeStress,
            const GeoType geoType
        );

    //- Destructor
    virtual ~neoHookeanElasticPlastic();


    // Member functions

        // Evaluation

            //- Pre update step
            virtual void preUpdate();

            //- Post update step
            virtual void postUpdate(const scalar f);

            //- Compute Cauchy stress and return wave speed
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
