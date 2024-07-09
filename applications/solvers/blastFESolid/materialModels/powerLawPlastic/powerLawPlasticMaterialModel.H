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
    Foam::materialModels::powerLawPlastic

Description
    Base class for a powerLawPlastic material

SourceFiles
    powerLawPlasticMaterialModel.C

\*---------------------------------------------------------------------------*/

#ifndef powerLawPlasticMaterialModel_H
#define powerLawPlasticMaterialModel_H

#include "elasticPlasticMaterialModel.H"
#include "plasticSolver.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace materialModels
{

/*---------------------------------------------------------------------------*\
                        Class powerLawPlastic Declaration
\*---------------------------------------------------------------------------*/

template<class Plastic>
class powerLawPlastic
:
    public Plastic
{
protected:
// Protected data

    //- Strength coefficient
    dimensionedScalar k_;

    //- Strain hardening exponent
    dimensionedScalar n_;

    //- Strain at yield
    dimensionedScalar epsilonY_;


    // Yield stress model
    struct YSModel
    :
        public plasticSolver::yieldStressModel
    {
        //- Strength coefficient
        const scalar& k_;

        //- Strain hardening exponent
        const scalar& n_;

        //- Strain at yield
        const scalar& epsilonY_;

    public:

        //- Constructor
        YSModel
        (
            const scalar& k,
            const scalar& n,
            const scalar& epsilonY
        )
        :
            k_(k),
            n_(n),
            epsilonY_(epsilonY)
        {}

        //- Return the current yield stress
        virtual scalar currentSigmaY
        (
            const scalar epsilonPEq0,
            const scalar curEpsilonPEq,
            const scalar J
        ) const
        {
            return k_*pow(curEpsilonPEq + epsilonY_, n_);
        }
    };


    //- Yield stress model
    YSModel ys_;


public:

    // Constructors

        //- Construct from order, size and element type
        powerLawPlastic
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
            k_("k", dimPressure, -1.0),
            n_("n", dimless, -1.0),
            epsilonY_("epsilonY", dimless, 0.0),
            ys_(k_.value(), n_.value(), epsilonY_.value())
        {
            scalar sigmaY
            (
                dict.lookupOrDefault<scalar>("sigmaY", -1.0)
            );
            scalar sigmaUT
            (
                dict.lookupOrDefault<scalar>("sigmaUT", -1.0)
            );
            scalar elongation
            (
                dict.lookupOrDefault<scalar>("elongation", -1.0)
            );

            k_.readIfPresent(dict);
            n_.readIfPresent(dict);

            if (n_.value() > 0)
            {}
            else if (sigmaY > 0 && sigmaUT > 0 && elongation > 0)
            {
                n_.value() = log(sigmaUT/sigmaY)/log(elongation/0.002);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Strain hardening coefficient could not be determined." << nl
                    << "Either provide n or sigmaY, sigmaUT, and elongation." << nl
                    << abort(FatalIOError);
            }

            if (k_.value() > 0)
            {}
            else if (sigmaUT > 0 && sigmaY > 0)
            {
                k_.value() = sigmaY/pow(sigmaY/this->E_.value() + 0.002, n_.value());
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Strength coefficient could not be determined." << nl
                    << "Either provide k or sigmaY and sigmaUT." << nl
                    << abort(FatalIOError);
            }

            if (sigmaY > 0)
            {
                epsilonY_.value() = pow(sigmaY/k_.value(), 1.0/n_.value());
            }
            else
            {
                epsilonY_.value() =
                    pow(this->E_.value()/k_.value(), 1.0/(n_.value() - 1.0));
            }

            if (!this->sigmaY_.headerOk())
            {
                this->sigmaY_ == k_*pow(epsilonY_, n_);
            }
        }


    //- Destructor
    virtual ~powerLawPlastic()
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
