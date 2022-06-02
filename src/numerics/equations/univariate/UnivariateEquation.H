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
    Foam::TypeEquation

Description
    Equation with multiple inputs and a single output

\*---------------------------------------------------------------------------*/

#ifndef UnivariateEquation_H
#define UnivariateEquation_H

#include "univariateEquation.H"
#include "multivariateEquation.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UnivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UnivariateEquation
:
    public univariateEquation<Type>
{
public:
    // Public typedef
    typedef univariateEquation<Type> BaseEquation;


protected:
// Protected data

    //- Lower limits
    scalarField lowerLimits_;

    //- Upper limit
    scalarField upperLimits_;

    //- Number in input variables
    const label nVar_;

    //- Finite difference step size
    mutable scalarField dX_;

    //- Calculate the gradient using a finite difference
    void calculateGradient
    (
        const scalarList& x0,
        const label li,
        List<Type>& grad
    ) const;


public:

    // Constructors

        //- Construct given number of equations
        UnivariateEquation(const label nVar, const string& eqnString = string::null);

        //- Construct given an equation string and bounds
        UnivariateEquation
        (
            const scalarList& lowerLimit,
            const scalarList& upperLimit,
            const string& eqnString = string::null
        );

        //- Construct from dictionary
        UnivariateEquation(const dictionary& dict);


    //- Destructor
    virtual ~UnivariateEquation();


    // Member Functions

        //- Return the number of input variables
        virtual label nVar() const
        {
            return nVar_;
        }

        //- Return the number of equations
        virtual label nEqns() const
        {
            return 1;
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
        virtual void setDX(const scalarList& newDx) const
        {
            dX_ = newDx;
        }

        //- Limit the value of x to the range of the equation
        virtual void limit(scalarField& x) const
        {
            x = max(min(upperLimits_, x), lowerLimits_);
        }


        //- Calculate the gradient
        virtual void dfdX
        (
            const scalarList& x,
            const label li,
            List<Type>& dfdx
        ) const;

        //- Calculate the gradient
        virtual void FX
        (
            const scalarList& x,
            const label li,
            List<Type>& fx
        ) const;

        //- Calculate the jacobian
        virtual void jacobian
        (
            const scalarList& x,
            const label li,
            List<Type>& fx,
            RectangularMatrix<Type>& J
        ) const;

        //- Check if the solution is bracketed
        virtual bool containsRoot
        (
            const List<Type>& y0s,
            const List<Type>& y1s
        ) const;

        //- Check if the solution is bracketed (compute bounds)
        virtual bool containsRoot(const label li) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UnivariateEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
