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
    Foam::univariateEquation

Description
    Abstract class for an equation that takes in a scalarList and outputs
    a single scalar value.

\*---------------------------------------------------------------------------*/

#ifndef univariateEquation_H
#define univariateEquation_H

#include "multivariateEquation.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class univariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class univariateEquation
:
    public multivariateEquation<Type>
{
public:

    // Constructors

        //- Null constructor
        univariateEquation()
        :
            multivariateEquation<Type>()
        {}

        //- Construct given name
        univariateEquation(const string& name)
        :
            multivariateEquation<Type>(name)
        {}


    //- Destructor
    virtual ~univariateEquation()
    {}

    //- Public member functions

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const = 0;

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const = 0;

        //- Return the finite difference step size
        virtual tmp<scalarField> dX() const = 0;

        //- Return the finite difference step size
        virtual void setDX(const scalarList& newDx) const = 0;

        //- Limit the value of x to the range of the equation
        virtual void limit(scalarField& x) const = 0;

        //- Return the function value
        virtual Type fX
        (
            const scalarList& x,
            const label li
        ) const = 0;

        //- Return the function value
        virtual void dfdX
        (
            const scalarList& x,
            const label li,
            List<Type>& dfdx
        ) const = 0;

        //- Return the function value
        virtual void d2fdX2
        (
            const scalarList& x,
            const label li,
            List<Type> d2dfx2
        ) const
        {
            NotImplemented;
        }

        //- Return the function value
        virtual void d3fdX3
        (
            const scalarList& x,
            const label li,
            List<Type> d3dfx3
        ) const
        {
            NotImplemented;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
