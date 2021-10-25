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
    Foam::Equation

Description
    Abstract base class for the an equation.

\*---------------------------------------------------------------------------*/

#ifndef Equation_H
#define Equation_H

#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
inline Type& setComponent(List<Type>& f, const direction diri)
{
    return f[diri];
}


template<class Type>
inline const Type& component(const List<Type>& f, const direction diri)
{
    return f[diri];
}

/*---------------------------------------------------------------------------*\
                          Class Equation Declaration
\*---------------------------------------------------------------------------*/

template<class InType, class OutType>
class Equation
{
protected:

    // Protected data

        //- Number in input variables
        const label nVar_;

        //- Number of equations
        const label nEqns_;

        //- Lower limits
        InType lowerLimits_;

        //- Upper limit
        InType upperLimits_;

        //- Finite difference step size
        InType dx_;


public:

    typedef Equation<InType, OutType> equationType;
    typedef InType inType;
    typedef OutType outType;

    // Constructors

        //- Construct given bounds
        Equation
        (
            const label nVar, 
            const label nEqns,
            const InType& lowerLimit, 
            const InType& upperLimit
        );


    //- Destructor
    virtual ~Equation();


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

        //- Return the number of derivatives
        virtual label nDerivatives() const = 0;

        //- Return the lower limit
        virtual const InType& lower() const
        {
            return lowerLimits_;
        }

        //- Access the lower limits
        virtual tmp<scalarField> lowerLimits() const;

        //- Access the upper limits
        virtual tmp<scalarField> upperLimits() const;

        //- Access the lower limit
        virtual InType& lower()
        {
            return lowerLimits_;
        }

        //- Return the upper limit
        virtual const InType& upper() const
        {
            return upperLimits_;
        }

        //- Access the upper limit
        virtual InType& upper()
        {
            return upperLimits_;
        }

        //- Return the finite difference step size for a direction
        virtual scalar dx(const label cmpti) const
        {
            return component(dx_, cmpti);
        }
        
        //- Return the finite difference step size for a direction
        virtual scalar& dx(const label cmpti)
        {
            return setComponent(dx_, cmpti);
        }

        //- Return the finite difference step size
        virtual tmp<scalarField> dx() const;

        //- Limit the value of x to the range of the equation
        virtual  void limit(InType& x) const
        {
            x = min(upperLimits_, max(x, lowerLimits_));
        }

        //- Return the function value
        virtual void f
        (
            const InType& x, 
            const label li,
            OutType& fx
        ) const 
        {
            NotImplemented;
        }
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
