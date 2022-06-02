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
    Foam::equation

Description
    Abstract base class an equation with a single input, and single output.
    Inherits from univariate and multivariate equations

\*---------------------------------------------------------------------------*/

#ifndef equation_H
#define equation_H

#include "univariateEquation.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class equation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class equation
:
    public univariateEquation<Type>
{
public:

    typedef scalar VarType;

    // Constructors

        //- Construct null
        equation()
        :
            univariateEquation<Type>()
        {}

        //- Construct from equation string
        equation(const string& eqnString)
        :
            univariateEquation<Type>(eqnString)
        {}

        //- Construct from dictionary
        equation(const dictionary& dict)
        :
            univariateEquation<Type>(dict)
        {}


    //- Destructor
    virtual ~equation()
    {}


    // Member Functions

        //- Return the lower limit
        virtual scalar lower() const = 0;

        //- Set the lower limit
        virtual void setLower(const scalar) = 0;

        //- Return the upper limit
        virtual scalar upper() const = 0;

        //- Set the upper limit
        virtual void setUpper(const scalar) = 0;

        //- Return the finite difference step size
        virtual scalar dx() const = 0;

        //- Set the finite difference step size
        virtual void setDx(const scalar& newDx) const = 0;

        //- Limit the value of x to the range of the equation
        using univariateEquation<Type>::limit;
        virtual void limit(scalar& x) const = 0;

        virtual Type fx(const scalar x, const label li) const = 0;

        //- Calculate the first derivative of the equation
        virtual Type dfdx(const scalar x, const label li) const
        {
            NotImplemented
            return this->fx(x, li);
        }

        //- Calculate the second derivative of the equation
        virtual Type d2fdx2(const scalar x, const label li) const
        {
            NotImplemented
            return this->fx(x, li);
        }

        //- Calculate the third derivative of the equation
        virtual Type d3fdx3(const scalar x, const label li) const
        {
            NotImplemented;
            return this->fx(x, li);
        }

        //- Calculate the fourth derivative of the equation
        virtual Type d4fdx4(const scalar x, const label li) const
        {
            NotImplemented;
            return this->fx(x, li);
        }

        //- Check if the solution is bracketed (compute bounds)
        using univariateEquation<Type>::containsRoot;
        virtual bool containsRoot(const Type& y0, const Type& y1) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
