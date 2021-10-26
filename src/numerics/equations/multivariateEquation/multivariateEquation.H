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
    Foam::scalarEquation

Description
    Abstract class for a multivariate equation

\*---------------------------------------------------------------------------*/

#ifndef multivariateEquation_H
#define multivariateEquation_H

#include "RectangularMatrix.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class multivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class multivariateEquation
{
public:

    // Constructors

        //- Construct given bounds
        multivariateEquation();


    //- Destructor
    virtual ~multivariateEquation();

    //- Public member functions

        //- Return the number of input variables
        virtual label nVar() const = 0;

        //- Return the number of equations
        virtual label nEqns() const = 0;

        //- Return the number of derivatives
        virtual label nDerivatives() const = 0;

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const = 0;

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const = 0;

        //- Return the finite difference step size
        virtual tmp<scalarField> dX() const = 0;

        //- Return the finite difference step size
        virtual void setDX(const scalarField& newDx) const = 0;

        //- Limit the value of x to the range of the equation
        virtual void limit(scalarField& x) const
        {
            x = max(lowerLimits(), min(x, upperLimits()));
        }

        //- Return the function value
        virtual void f
        (
            const scalarField& x,
            const label li,
            Field<Type>& fx
        ) const = 0;

        //- Return the function value
        virtual void jacobian
        (
            const scalarField& x,
            const label li,
            Field<Type>& fx,
            RectangularMatrix<Type>& J
        ) const = 0;

        //- Check if the solution is bracketed
        virtual bool containsRoot
        (
            const Field<Type>& y0s,
            const Field<Type>& y1s
        ) const;

        //- Check if the solution is bracketed (compute bounds)
        virtual bool containsRoot(const label li) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "multivariateEquation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
