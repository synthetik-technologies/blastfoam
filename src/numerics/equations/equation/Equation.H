/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::Equation

Description
    Equation with a single input and single output. Derived from equation,
    univariateEquation and multivariateEquation so that an equation can be
    cast in any form and used in any of the multivariate numerical methods.

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
public:
    // Public typedef
    typedef equation<Type> BaseEquation;


protected:
    // Protected Data

        //- Lower limit
        scalar lower_;

        //- Upper limit
        scalar upper_;

        //- Finite difference step size
        mutable scalar dx_;


public:

    // Constructors

        //- Construct null
        Equation(const string& eqnString = string::null);

        //- Construct given bounds and name
        Equation
        (
            const scalar lower,
            const scalar upper,
            const string& eqnString = string::null
        );

         //- Construct from dictionary
        Equation(const dictionary& dict);


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

        //- Set the lower limit
        virtual void setLower(const scalar low)
        {
            lower_ = low;
        }

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const
        {
            return tmp<scalarField>(new scalarField(1, lower_));
        }

        //- Set the lower limit
        virtual void setLowerLimits(const scalarList& low)
        {
            lower_ = low[0];
        }

        //- Return the upper limit
        virtual scalar upper() const
        {
            return upper_;
        }

        //- Set the upper limit
        virtual void setUpper(const scalar up)
        {
            upper_ = up;
        }

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const
        {
            return tmp<scalarField>(new scalarField(1, upper_));
        }

        //- Set the lower limit
        virtual void setUpperLimits(const scalarList& up)
        {
            upper_ = up[0];
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
            const UList<scalar>& x,
            const label li
        ) const;

        //- Calulate the function at a point
        virtual void FX
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& fx
        ) const;

        //- Calulate the gradient
        virtual void dfdX
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& dfdx
        ) const;

        //- Calulate the Jacobian (dfdx)
        virtual void jacobian
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
        ) const;

        //- Check if the a root is on the interval (change in signs)
        virtual bool containsRoot(const label li) const;
        virtual bool containsRoot(const Type& y0, const Type& y1) const;
        virtual bool containsRoot
        (
            const UList<Type>& y0s,
            const UList<Type>& y1s
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
