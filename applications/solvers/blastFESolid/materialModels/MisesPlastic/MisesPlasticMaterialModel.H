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
    Foam::materialModels::MisesPlastic

Description
    Base class for a MisesPlastic material

SourceFiles
    MisesPlasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef MisesPlasticMaterialModel_H
#define MisesPlasticMaterialModel_H

#include "elasticPlasticMaterialModel.H"
#include "plasticSolver.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class MisesPlastic Declaration
\*---------------------------------------------------------------------------*/

template<class Plastic>
class MisesPlastic
:
    public Plastic
{
protected:
// Protected data

    //- An iterative procedure is used when the yield stress is a nonlinear
    //  function of plastic strain
    const bool nonLinearPlasticity_;

    //- Table of post-yield stress versus plastic strain
    autoPtr<Function1<scalar>> stressPlasticStrainSeries_;

    //- Linear plastic modulus. It is only used if plasticity is linear,
    //  defined by two points on the stress-plastic strain curve
    scalar Hp_;

    // Yield stress model
    struct YSModel
    :
        public plasticSolver::yieldStressModel
    {
    private:

        // Pointer to stess-plastic strain function
        const Function1<scalar>* stressPlasticStrainSeries_;

    public:
        YSModel(const autoPtr<Function1<scalar>>& spss)
        :
            stressPlasticStrainSeries_
            (
                spss.valid()
              ? &spss()
              : nullptr
            )
        {}

        //- Return the current yield stress
        virtual scalar currentSigmaY
        (
            const scalar epsilonPEq0,
            const scalar curEpsilonPEq,
            const scalar J
        ) const
        {
            return
                stressPlasticStrainSeries_->value
                (
                    max(epsilonPEq0, 0.0)
                )*J;
        }
    };


    //- Yield stress model
    YSModel ys_;


public:

    // Constructors

        //- Construct from order, size and element type
        MisesPlastic
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
            nonLinearPlasticity_
            (
                dict.lookupOrDefault("nonLinearPlasticity", false)
            ),
            stressPlasticStrainSeries_
            (
                nonLinearPlasticity_
              ? Function1<scalar>::New
                (
                    "stressPlasticStrainSeries",
                    mesh.time().userUnits(),
                    dimPressure,
                    dict
                ).ptr()
              : nullptr
            ),
            Hp_
            (
                nonLinearPlasticity_
              ? 0.0
              : dict.lookup<scalar>("Hp")
            ),
            ys_(stressPlasticStrainSeries_)
        {
            // if (!this->sigmaY_.headerOk())
            {
                dimensionedScalar sigmaY("sigmaY", dimPressure, 0.0);
                if (nonLinearPlasticity_)
                {
                    sigmaY.value() = stressPlasticStrainSeries_->value(0.0);
                }
                else
                {
                    sigmaY.read(dict);
                }

                this->sigmaY_ = sigmaY.value();
            }
        }


    //- Destructor
    virtual ~MisesPlastic()
    {}


    // Member functions

        // Evaluation

            //- Return a yield stress model
            virtual tmp<plasticSolver::yieldStressModel> yieldStress
            (
                const label elemi,
                const label rulei
            ) const
            {
                return ys_;
            }


    // Member functions

        //- Update plasticN, DLambda, DSigmaY and sigmaY
        virtual void updatePlasticity
        (
            const plasticSolver::yieldStressModel& ys,
            symmTensor& plasticN,
            scalar& DLambda,
            scalar& sigmaY,
            const scalar sigmaY0,
            const scalar fTrial,
            const symmTensor& sTrial,
            const scalar epsilonPEq0,
            const scalar muBar,
            const scalar maxMagDEpsilon,
            const scalar J
        ) const
        {
            if (fTrial < small)
            {
                plasticN = symmTensor::I;
                DLambda = 0.0;
                return;
            }

            const scalar magS = mag(sTrial);
            if (magS > small)
            {
                plasticN = sTrial/magS;
            }
            else
            {
                plasticN = symmTensor::I;
            }

            if (!nonLinearPlasticity_)
            {
                DLambda = fTrial/(2.0*muBar);
                if (mag(Hp_) > small)
                {
                    DLambda /= (1.0 + Hp_/(3.0*muBar));
                    sigmaY = sigmaY0 + DLambda*Hp_;
                }
                return;
            }
            this->newtonLoop
            (
                ys,
                DLambda,
                sigmaY,
                epsilonPEq0,
                magS,
                muBar,
                maxMagDEpsilon,
                J
            );
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
