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
    Foam::materialModels::JohnsonCookPlastic

Description
    Base class for a JohnsonCookPlastic material

SourceFiles
    JohnsonCookPlasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef JohnsonCookPlasticMaterialModel_H
#define JohnsonCookPlasticMaterialModel_H

#include "elasticPlasticMaterialModel.H"
#include "plasticSolver.H"
#include "volPointInterpolation.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class JohnsonCookPlastic Declaration
\*---------------------------------------------------------------------------*/

template<class Plastic>
class JohnsonCookPlastic
:
    public Plastic
{
protected:
// Protected data

    //- Yield stress at reference conditions
        dimensionedScalar A_;

        //- Strain hardening constant
        dimensionedScalar B_;

        //- Strengthing coefficient of strain rate
        dimensionedScalar C_;

        //- Strain hardening exponentn
        dimensionedScalar n_;

        //- Reference plastic strain
        dimensionedScalar epsilonP0_;

        //- Referenece plastic strain rate
        dimensionedScalar epsilonPDot0_;

        //- Is the total strain rate used
        Switch totalStrainRate_;

        //- Are temperature effects included
        Switch temperatureEffects_;

        //- Thermal softening coefficient
        dimensionedScalar m_;

        //- Reference temperature
        dimensionedScalar Tref_;

        //- Melting temperature
        dimensionedScalar Tm_;


    // Yield stress model
    struct YSModel
    :
        public plasticSolver::yieldStressModel
    {
        const scalar& A_;
        const scalar& B_;
        const scalar& C_;
        const scalar& n_;
        const scalar& epsilonP0_;
        const scalar& epsilonPDot0_;
        const bool& totalStrainRate_;
        const bool& temperatureEffects_;
        const scalar& m_;
        const scalar& Tref_;
        const scalar& Tm_;

        const scalar dt_;
        const scalar T_;
        const scalar epsilonDot_;

    public:

        //- Constructor
        YSModel
        (
            const scalar& A,
            const scalar& B,
            const scalar& C,
            const scalar& n,
            const scalar& epsilonP0,
            const scalar& epsilonPDot0,
            const bool& totalStrainRate,
            const bool& temperatureEffects,
            const scalar& m,
            const scalar& Tref,
            const scalar& Tm
        )
        :
            A_(A),
            B_(B),
            C_(C),
            n_(n),
            epsilonP0_(epsilonP0),
            epsilonPDot0_(epsilonPDot0),
            totalStrainRate_(totalStrainRate),
            temperatureEffects_(temperatureEffects),
            m_(m),
            Tref_(Tref),
            Tm_(Tm),
            dt_(great),
            T_(Tref),
            epsilonDot_(0.0)
        {}
        YSModel
        (
            const YSModel& ys,
            const scalar dt,
            const scalar T,
            const scalar epsilonDot
        )
        :
            A_(ys.A_),
            B_(ys.B_),
            C_(ys.C_),
            n_(ys.n_),
            epsilonP0_(ys.epsilonP0_),
            epsilonPDot0_(ys.epsilonPDot0_),
            totalStrainRate_(ys.totalStrainRate_),
            temperatureEffects_(ys.temperatureEffects_),
            m_(ys.m_),
            Tref_(ys.Tref_),
            Tm_(ys.Tm_),
            dt_(dt),
            T_(T),
            epsilonDot_(epsilonDot)
        {}

        //- Return the current yield stress
        virtual scalar currentSigmaY
        (
            const scalar epsilonPEq0,
            const scalar curEpsilonPEq,
            const scalar J
        ) const
        {
            scalar epsilonDot
            (
                (curEpsilonPEq - epsilonPEq0)/dt_
            );
            if (totalStrainRate_)
            {
                epsilonDot += epsilonDot_;
            }
            scalar sigmay =
                (A_ + B_*pow(max(curEpsilonPEq/epsilonP0_, 0), n_))
               *(
                    1.0
                  + (C_*log(max(epsilonDot/epsilonPDot0_, 1.0)))
                );
            if (temperatureEffects_)
            {
                sigmay *= (1.0 - (T_ - Tref_)/(Tm_ - Tref_));
            }
            return J*sigmay;
        }
    };

    //- Yield stress model
    YSModel ys_;

    //- Temperature field
    tmp<pointScalarField> TPtr_;


public:

    // Constructors

        //- Construct from order, size and element type
        JohnsonCookPlastic
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
            A_("A", dimPressure, dict),
            B_("B", dimPressure, dict),
            C_("C", dimless, dict),
            n_("n", dimless, dict),
            epsilonP0_
            (
                dimensionedScalar::lookupOrDefault
                (
                    "epsilonP0",
                    dict,
                    dimless,
                    1.0
                )
            ),
            epsilonPDot0_
            (
                dimensionedScalar::lookupOrDefault
                (
                    "epsilonPDot0",
                    dict,
                    inv(dimTime),
                    1.0
                )
            ),
            totalStrainRate_
            (
                dict.lookupOrDefault<Switch>("totalStrainRate", false)
            ),
            temperatureEffects_
            (
                dict.lookupOrDefault<Switch>("temperatureEffects", false)
            ),
            m_(1.0),
            Tref_("Tref", dimTemperature, 300),
            Tm_("Tm", dimTemperature, 1000),
            ys_
            (
                A_.value(),
                B_.value(),
                C_.value(),
                n_.value(),
                epsilonP0_.value(),
                epsilonPDot0_.value(),
                totalStrainRate_,
                temperatureEffects_,
                m_.value(),
                Tref_.value(),
                Tm_.value()
            )
        {
            if (temperatureEffects_)
            {
                m_ = dict.lookup<scalar>("m");
                Tref_.read(dict);
                Tm_.read(dict);
            }

            if (!this->sigmaY_.headerOk())
            {
                this->sigmaY_ == A_;
            }
        }


    //- Destructor
    virtual ~JohnsonCookPlastic()
    {}


    // Member functions

        //- Return the yield stress model
        virtual tmp<plasticSolver::yieldStressModel> yieldStress
        (
            const label elemi,
            const label rulei
        ) const
        {
            if (!totalStrainRate_ && !temperatureEffects_)
            {
                return ys_;
            }
            const element& elem = this->mesh_.elements()[elemi];

            scalar epsilonDot = 0.0;
            scalar T = Tref_.value();
            if (totalStrainRate_)
            {
                epsilonDot =
                    mag
                    (
                        dev
                        (
                            symm
                            (
                                this->gradient
                                (
                                    UIndirectList<vector>(this->U_, elem),
                                    this->mesh_.dshapes()[elemi][rulei],
                                    this->mesh_.invJs()[elemi][rulei]
                                )
                            )
                        )
                    )*sqrt(2.0/3.0);
            }
            if (temperatureEffects_)
            {
                UIndirectList<scalar> T_loc(TPtr_(), elem);
                const scalarList& shape = this->mesh_.shapes()[elemi][rulei];
                forAll(shape, si)
                {
                    T += T_loc[si]*shape[si];
                }
            }

            return tmp<plasticSolver::yieldStressModel>
            (
                new YSModel
                (
                    ys_,
                    this->mesh_.time().deltaTValue(),
                    T,
                    epsilonDot
                )
            );
        }

        //- Pre update step
        virtual void preUpdate()
        {
//             if (temperatureEffects_)
//             {
//                 const polyMesh& mesh = this->mesh_.mesh();
//                 if (mesh.foundObject<pointScalarField>("pointT"))
//                 {
//                     if (TPtr_.empty())
//                     {
//                         TPtr_ = tmp<pointScalarField>
//                         (
//                             mesh.lookupObject<pointScalarField>("pointT")
//                         );
//                     }
//                 }
//                 else
//                 {
//                     TPtr_ =
//                         volPointInterpolation::New(mesh).interpolate
//                         (
//                             mesh.template lookupObject<volScalarField>("T")
//                         );
//                 }
//             }
            Plastic::preUpdate();
        }

        //- Post update step
        virtual void postUpdate(const scalar f)
        {
            if (TPtr_.isTmp())
            {
                TPtr_.clear();
            }
            Plastic::postUpdate(f);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace materialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
