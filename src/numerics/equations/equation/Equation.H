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
    Foam::equation

Description
    Scalar equation with a single input. Derived from both scalarEquation
    and multivariateEquation so that a scalar equation can be cast in either
    form.

\*---------------------------------------------------------------------------*/

#ifndef Equation_H
#define Equation_H

#include "equation.H"
#include "univariateEquation.H"
#include "multivariateEquation.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class Equation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Equation
:
    public equation<Type>
{
protected:
    // Protected Data

        //- Lower limits
        scalar lower_;

        //- Upper limit
        scalar upper_;

        //- Finite difference step size
        mutable scalar dx_;


public:

    // Constructors

        //- Construct given bounds
        Equation(const scalar lower, const scalar upper);

        //- Construct given bounds and name
        Equation(const string& name, const scalar lower, const scalar upper);


    //- Destructor
    virtual ~Equation();


    // Member Functions

        //- Return the number of input variables
        virtual label nVar() const
        {
            return 1;
        }

        //- Return the number of equations
        virtual label nEqns() const
        {
            return 1;
        }

        //- Return the lower limit
        virtual scalar lower() const
        {
            return lower_;
        }

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const
        {
            return tmp<scalarField>(new scalarField(1, lower_));
        }

        //- Return the upper limit
        virtual scalar upper() const
        {
            return upper_;
        }

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const
        {
            return tmp<scalarField>(new scalarField(1, upper_));
        }

        //- Return the finite difference step size
        virtual scalar dx() const
        {
            return dx_;
        }

        //- Return the finite difference step size
        virtual tmp<scalarField> dX() const
        {
            return tmp<scalarField>(new scalarField(1, dx_));
        }

        //- Return the finite difference step size
        virtual void setDx(const scalar& newDx) const
        {
            dx_ = newDx;
        }

        //- Return the finite difference step size
        virtual void setDX(const scalarList& newDx) const
        {
            dx_ = newDx[0];
        }

        //- Limit the value of x to the range of the equation
        virtual void limit(scalar& x) const
        {
            x = max(min(upper_, x), lower_);
        }
        virtual void limit(scalarField& x) const
        {
            limit(x[0]);
        }

        //- Calulate the function at a point
        virtual Type fX
        (
            const scalarList& x,
            const label li
        ) const;

        //- Calulate the function at a set of points
        virtual void FX
        (
            const scalarList& x,
            const label li,
            List<Type>& fx
        ) const;

        //- Calulate the gradient
        virtual tmp<Field<Type>> dfdX
        (
            const scalarList& x,
            const label li
        ) const;

        //- Calulate the Jacobian (dfdx)
        virtual void jacobian
        (
            const scalarList& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
        ) const;

        virtual bool containsRoot(const label li) const;
        virtual bool containsRoot(const Type& y0, const Type& y1) const;
        virtual bool containsRoot
        (
            const List<Type>& y0s,
            const List<Type>& y1s
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Equation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
