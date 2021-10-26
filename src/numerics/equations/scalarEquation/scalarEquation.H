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
    Base class for the a scalar equation

\*---------------------------------------------------------------------------*/

#ifndef scalarEquation_H
#define scalarEquation_H

#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class scalarEquation Declaration
\*---------------------------------------------------------------------------*/

class scalarEquation
{
public:

    // Constructors

        //- Construct given bounds
        scalarEquation()
        {}


    //- Destructor
    virtual ~scalarEquation()
    {}

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
            scalar& fx
        ) const = 0;

        //- Return the function value
        virtual void gradient
        (
            const scalarField& x,
            const label li,
            scalar& fx,
            scalarField& grad
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
