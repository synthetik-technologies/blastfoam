/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
    Foam::MultivariateEquation

Description
    Abstract base class a multivariate equation

\*---------------------------------------------------------------------------*/

#ifndef MultivariateEquation_H
#define MultivariateEquation_H

#include "scalarField.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MultivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MultivariateEquation
{
protected:

    // Protected data

        //- Lower limits
        scalarField lowerLimits_;

        //- Upper limits
        scalarField upperLimits_;

        //- Finite difference step size
        scalarField dx_;


    // Protected member functions

        //- Check if derivatives have been implemented
//         bool checkJacobian() const;
//         bool checkHessian() const;

        //- Make sure limits have been set
        void checkLimits() const;

        //- Create a Jacobian matrix given a list of dx
        void calculateJacobian
        (
            const scalarList& x0,
            const label li,
            const List<Type>& f0,
            RectangularMatrix<Type>& J
        ) const;

public:

    // Constructors

        //- Construct null
        MultivariateEquation(const label nEqns);

        //- Construct with limits
        MultivariateEquation
        (
            const scalarList& lowerLimits,
            const scalarList& upperLimits
        );


    //- Destructor
    virtual ~MultivariateEquation();


    // Member Functions

        //- Return the number of equations in the system
        virtual label nEqns() const = 0;

        //- Return the number of derivatives in the system
        virtual label nDerivatives() const = 0;

        //- Return the lower limits
        const scalarList& lower() const
        {
            return lowerLimits_;
        }

        //- Access the lower limits
        scalarList& lowerLimits()
        {
            return lowerLimits_;
        }

        //- Return the upper limits
        const scalarList& upper() const
        {
            return upperLimits_;
        }

        //- Access the upper limits
        scalarList& upper()
        {
            return upperLimits_;
        }

        //- Return the finite difference step size
        const scalarList& dx() const
        {
            return dx_;
        }

        //- Access the finite difference step size
        scalarList& dx()
        {
            return dx_;
        }

        //- Limit the value of x to the range of the equation
        inline void limit(scalarList& x) const
        {
            x = min(upperLimits_, max(x, lowerLimits_));
        }

        //- Check if the coordinates are with in the range
        bool checkBounds(const scalarList& xs) const;

        //- Return the function value
        virtual void f
        (
            const scalarList& x,
            const label li,
            List<Type>& fx
        ) const = 0;

        //- Return the function value
        tmp<Field<Type>> f
        (
            const scalarList& x,
            const label li
        ) const;

        //- Calculate the Jacobian of the equation
        virtual void jacobian
        (
            const scalarList& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& dfdx
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "MultivariateEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
