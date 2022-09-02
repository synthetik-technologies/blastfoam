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
    Foam::MultivariateEquation

Description
    Equation for multiple inputs and outputs.

\*---------------------------------------------------------------------------*/

#ifndef MultivariateEquation_H
#define MultivariateEquation_H


#include "equationBase.H"
#include "scalarField.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class multivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateEquation
:
    public equationBase
{
public:

    typedef UList<scalar> VarType;
    typedef List<Type> EqnType;
    typedef Type CmptType;

    //- Constructors

        //- Null constructor
        multivariateEquation()
        :
            equationBase()
        {}

        //- Construct from equation string
        multivariateEquation(const string& eqnString)
        :
            equationBase(eqnString)
        {}

        //- Construct from a list of equation strings
        multivariateEquation(const List<string>& eqnStrings)
        :
            equationBase(eqnStrings)
        {}

        //- Construct from dictionary
        multivariateEquation(const dictionary& dict)
        :
            equationBase(dict)
        {}


    //- Destructor
    virtual ~multivariateEquation()
    {}


    //- Public member functions

        //- Return the number of input variables
        virtual label nVar() const = 0;

        //- Return the number of equations
        virtual label nEqns() const = 0;

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const = 0;

        //- Set the lower limit
        virtual void setLowerLimits(const scalarList&) = 0;

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const = 0;

        //- Set the upper limit
        virtual void setUpperLimits(const scalarList&) = 0;

        //- Return the finite difference step size
        virtual tmp<scalarField> dX() const = 0;

        //- Return the finite difference step size
        virtual void setDX(const scalarList& newDx) const = 0;

        //- Limit the value of x to the range of the equation
        virtual void limit(scalarField& x) const
        {
            x = max(lowerLimits(), min(x, upperLimits()));
        }

        //- Return the function value
        virtual void FX
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& fx
        ) const = 0;

        //- Return the function value
        virtual void jacobian
        (
            const UList<scalar>& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
        ) const = 0;

        //- Check if the solution is bracketed
        virtual bool containsRoot
        (
            const UList<Type>& y0s,
            const UList<Type>& y1s
        ) const = 0;

        //- Check if the solution is bracketed (compute bounds)
        virtual bool containsRoot(const label li) const = 0;
};


/*---------------------------------------------------------------------------*\
                          Class MultivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MultivariateEquation
:
    public multivariateEquation<Type>
{
public:
    // Public typedef
    typedef multivariateEquation<Type> BaseEquation;

protected:
    // Protected data

        //- Lower limits
        scalarField lowerLimits_;

        //- Upper limit
        scalarField upperLimits_;

        //- Number in input variables
        const label nVar_;

        //- Number of equations
        const label nEqns_;

        //- Finite difference step size
        mutable scalarField dX_;


    // Protected member functions

        //- Create a Jacobian matrix given a list of dx
        void calculateJacobian
        (
            const UList<scalar>& x0,
            const label li,
            const List<Type>& f0,
            RectangularMatrix<Type>& J
        ) const;


public:

    // Constructors

        //- Construct with number of equations and number of variabels
        MultivariateEquation
        (
            const label nEqns,
            const label nVar,
            const List<string>& eqnString = List<string>()
        );

        //- Construct with number of equations, limits, and name
        MultivariateEquation
        (
            const label nEqns,
            const scalarList& lowerLimits,
            const scalarList& upperLimits,
            const List<string>& eqnString = List<string>()
        );

        //- Construct from dictionary
        MultivariateEquation(const dictionary& dict);


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

        //- Set the lower limit
        virtual void setLowerLimits(const scalarList& low)
        {
            lowerLimits_ = low;
        }

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const
        {
            return upperLimits_;
        }

        //- Set the upper limit
        virtual void setUpperLimits(const scalarList& up)
        {
            upperLimits_ = up;
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
            const UList<scalar>& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
        ) const;

        //- Check if the solution is bracketed
        virtual bool containsRoot
        (
            const UList<Type>& y0s,
            const UList<Type>& y1s
        ) const;

        //- Check if the solution is bracketed (compute bounds)
        virtual bool containsRoot(const label li) const;
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
