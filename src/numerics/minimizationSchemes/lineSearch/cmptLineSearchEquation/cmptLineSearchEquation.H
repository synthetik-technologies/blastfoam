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
    Foam::cmptLineSearchEquation

Description
    Equation for solving a multivariate, single output equation given a
    reference point, direction, and the distance as an input

SourceFiles
    cmptLineSearchEquation.C

\*---------------------------------------------------------------------------*/

#ifndef cmptLineSearchEquation_H
#define cmptLineSearchEquation_H

#include "UnivariateEquationsFwd.H"
#include "EquationsFwd.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class univariateMinimizationScheme;

/*---------------------------------------------------------------------------*\
                Class cmptLineSearchEquation Declaration
\*---------------------------------------------------------------------------*/

class cmptLineSearchEquation
:
    public ScalarEquation
{
private:
    // Private member data

        //- Reference to the actual equation
        const scalarUnivariateEquation& eqns_;

        //- Scalar Equation (optional)
        UautoPtr<const scalarEquation> sEqn_;

        //- Travel component
        label cmpt_;

        //- Sign of travel direction
        label sign_;

        //- Reference values
        scalarField x0_;

        //- tmp class for storage of the traversal x field
        mutable scalarField x_;

        //- Tmp field for dfdX
        mutable scalarField dfdX_;

        //- Line search minimization scheme
        mutable autoPtr<univariateMinimizationScheme> lineSearcher_;

        // Bracketing controls

            //- Maximum number of bracketing iterations
            label maxBracketIter_;

            //- Initial distance
            scalar dx0_;

            //- Step size increase factor
            scalar k_;


public:

    // Constructor

        // Construct from the given equation
        cmptLineSearchEquation
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        // Construct from the given equation
        cmptLineSearchEquation
        (
            const scalarUnivariateEquation& eqn,
            const cmptLineSearchEquation& ls
        );

    // Destructor
    virtual ~cmptLineSearchEquation();


    // Member Functions

        // Access

            //- Return the number of derivatives
            inline label nDerivatives() const
            {
                return eqns_.nDerivatives();
            }


        // Edit

            //- Set direction based on gradient and start of line search
            void setX0(const scalarList& x0);

            //- Set travel component
            void setCmpt(const label cmpt, const label sign = 1);

            //- Flip travel direction
            void flip();


            //- Return the true x values given a distance from the
            //  reference state
            const scalarField& calcX(const scalar dist) const;


        // Evaluation

            //- Conduct line search
            void search
            (
                const label li,
                scalarList& xNew
            );


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
                return dfdX_[cmpt_];
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
