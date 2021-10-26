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
        mutable InType dx_;


public:

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
        inline label nVar() const;

        //- Return the number of equations
        inline label nEqns() const;

        //- Return the number of derivatives
        virtual label nDerivatives() const = 0;

        //- Return the lower limits
        virtual const InType& lower() const;

        //- Access the lower limits
        virtual InType& lower();

        //- Return the upper limits
        inline const InType& upper() const;

        //- Access the upper limits
        inline InType& upper();

        //- Return the finite difference step size
        inline const InType& dx() const;

        //- Access the finite difference step size
        inline InType& dx();

        //- Return the finite difference step size for a direction
        inline scalar dx(const label cmpti) const;

        //- Limit the value of x to the range of the equation
        inline  void limit(InType& x) const;

        //- Return the function value
        virtual void f
        (
            const InType& x,
            const label li,
            OutType& fx
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EquationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "Equation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
