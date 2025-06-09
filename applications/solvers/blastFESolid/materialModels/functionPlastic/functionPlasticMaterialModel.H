/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
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
    Foam::materialModels::functionPlastic

Description
    Base class for a functionPlastic material

SourceFiles
    functionPlasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef functionPlasticMaterialModel_H
#define functionPlasticMaterialModel_H

#include "elasticPlasticMaterialModel.H"
#include "plasticSolver.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class functionPlastic Declaration
\*---------------------------------------------------------------------------*/

template<class Plastic>
class functionPlastic
:
    public Plastic
{
protected:
// Protected data

    // Yield stress function
    autoPtr<Function1<scalar>> stressStrain_;


    // Yield stress model
    struct YSModel
    :
        public plasticSolver::yieldStressModel
    {
        //- Yield stress function
        const Function1<scalar>& stressStrain_;

    public:

        //- Constructor
        YSModel(const Function1<scalar>& stressStrain)
        :
            stressStrain_(stressStrain)
        {}

        //- Return the current yield stress
        virtual scalar currentSigmaY
        (
            const scalar epsilonPEq0,
            const scalar curEpsilonPEq,
            const scalar J
        ) const
        {
            return J*stressStrain_.value(max(curEpsilonPEq, small));
        }
    };


    //- Yield stress model
    YSModel ys_;


public:

    // Constructors

        //- Construct from order, size and element type
        functionPlastic
        (
            const dictionary& dict,
            const feMesh1& mesh,
            const pointVectorField& D,
            const pointVectorField& U,
            const bool planeStress,
            const materialModel::GeoType geoType
        )
        :
            Plastic(dict, mesh, D, U, planeStress, geoType),
            stressStrain_
            (
                Function1<scalar>::New
                (
                    "stressPlasticStrainSeries",
                    mesh.time().userUnits(),
                    dimPressure,
                    dict
                )
            ),
            ys_(stressStrain_())
        {
            // if (!this->sigmaY_.headerOk())
            {
                this->sigmaY_ =
                    dimensionedScalar
                    (
                        "sigmaY",
                        dimPressure,
                        stressStrain_->value(0.0)
                    ).value();
            }
        }


    //- Destructor
    virtual ~functionPlastic()
    {}


    // Member functions

        //- Return the yield stress model
        virtual tmp<plasticSolver::yieldStressModel> yieldStress
        (
            const label elemi,
            const label rulei
        ) const
        {
            return ys_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
