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
    Abstract base class a multivariate equation. Uses a scalar field as
    the input, but can output a list of any type.

\*---------------------------------------------------------------------------*/

#ifndef MultivariateEquation_H
#define MultivariateEquation_H

#include "multivariateEquation.H"
#include "scalarField.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MultivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MultivariateEquation
:
    public multivariateEquation<Type>
{
protected:
    // Protected data

        //- Number in input variables
        const label nVar_;

        //- Number of equations
        const label nEqns_;

        //- Lower limits
        scalarField lowerLimits_;

        //- Upper limit
        scalarField upperLimits_;

        //- Finite difference step size
        mutable scalarField dX_;


    // Protected member functions

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

        //- Construct with number of equations and limits
        MultivariateEquation
        (
            const label nEqns,
            const scalarList& lowerLimits,
            const scalarList& upperLimits
        );

        //- Construct with number of equations, limits, and name
        MultivariateEquation
        (
            const string& name,
            const label nEqns,
            const scalarList& lowerLimits,
            const scalarList& upperLimits
        );


    //- Destructor
    virtual ~MultivariateEquation();


    // Member Functions

        //- Return the number of input variables
        virtual label nVar() const
        {
            return nVar_;
        }

        //- Return the number of equations
        virtual label nEqns() const
        {
            return nEqns_;
        }

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const
        {
            return lowerLimits_;
        }

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const
        {
            return upperLimits_;
        }

        //- Return the finite difference step size
        virtual tmp<scalarField> dX() const
        {
            return dX_;
        }

        //- Return the finite difference step size
        virtual void setDX(const scalarList& newDX) const
        {
            dX_ = newDX;
        }

        //- Limit the value of x to the range of the equation
        virtual void limit(scalarField& x) const
        {
            x = max(min(upperLimits_, x), lowerLimits_);
        }

        //- Calculate the Jacobian of the equation
        virtual void jacobian
        (
            const scalarList& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
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
