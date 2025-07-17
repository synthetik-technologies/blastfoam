/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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
    Foam::lineSearchEquation

Description
    Equation for solving a multivariate, single output equation given a
    reference point, direction, and the distance as an input

SourceFiles
    lineSearchEquation.C

\*---------------------------------------------------------------------------*/

#ifndef lineSearchEquation_H
#define lineSearchEquation_H

#include "UnivariateEquationsFwd.H"
#include "EquationsFwd.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class lineSearchEquation Declaration
\*---------------------------------------------------------------------------*/

class lineSearchEquation
:
    public ScalarEquation
{
public:
    // Public enumeration

        enum DescentMethod
        {
            NONE,
            FLETCHER_REEVES,
            POLAK_RIBIERE
        };

        static const NamedEnum<DescentMethod, 3> DescentMethodNames_;


private:
    // Private member data

        //- Reference to the actual equation
        const scalarUnivariateEquation& eqns_;

        //- Scalar Equation (optional)
        UautoPtr<const scalarEquation> sEqn_;

        //- Has this equation been initialized
        bool initialized_;

        //- Descent type
        DescentMethod descent_;

        //- Copy of the gradient
        scalarField grad_;

        //- Travel direction
        scalarField dir_;

        //- Reference values
        scalarField x0_;

        //- Descent coefficient
        scalar beta_;

        //- tmp class for storage of the traversal x field
        mutable scalarField x_;

        //- Tmp field for dfdX
        mutable scalarField dfdX_;


public:

    // Constructor

        // Construct from the given equation
        lineSearchEquation
        (
            const scalarUnivariateEquation& eqn,
            const word& descentType
        );

        // Construct from the given equation and descent type
        lineSearchEquation
        (
            const scalarUnivariateEquation& eqn,
            const DescentMethod descent
        );

        // Construct from the given equation and descent type
        lineSearchEquation
        (
            const scalarUnivariateEquation& eqn,
            const lineSearchEquation& ls
        );


    // Destructor
    virtual ~lineSearchEquation();


    // Member Functions

        // Access

            //- Return the number of derivatives
            inline label nDerivatives() const
            {
                return eqns_.nDerivatives();
            }

            //- Return the travel direction
            inline const scalarList& dir() const
            {
                return dir_;
            }

            //- Return the descent coefficient
            inline scalar beta() const
            {
                return beta_;
            }


        // Edit

            //- Set direction based on gradient and start of line search
            void update(const scalarList& x0, const scalarList& grad);

            //- Set direction based on gradient and start of line search
            void updateDir(const scalarList& x0, const scalarList& dir);

            //- Return the true x values given a distance from the
            //  reference state
            const scalarField& calcX(const scalar dist) const;


        // Evaluation

            //- Return the value of the function by advancing the given distance
            //  from the reference point
            inline scalar fx(const scalar dist, const label li) const
            {
                return eqns_.fX(calcX(dist), li);
            }

            //- Return the first derivative
            inline scalar dfdx(const scalar dist, const label li) const
            {
                if (sEqn_.valid())
                {
                    return sEqn_->dfdx(calcX(dist)[0], li);
                }
                eqns_.dfdX(calcX(dist), li, dfdX_);

                scalar dfdx = 0.0;
                forAll(dir_, i)
                {
                    dfdx += dfdX_[i]*dir_[i];
                }
                return dfdx;
            }

            //- Return the second derivative
            inline scalar d2fdx2(const scalar dist, const label li) const
            {
                if (sEqn_.valid())
                {
                    return sEqn_->d2fdx2(calcX(dist)[0], li);
                }
                NotImplemented;
                return 0;
            }

            //- Return the third derivative
            inline scalar d3fdx3(const scalar dist, const label li) const
            {
                if (sEqn_.valid())
                {
                    return sEqn_->d3fdx3(calcX(dist)[0], li);
                }
                NotImplemented;
                return 0;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
