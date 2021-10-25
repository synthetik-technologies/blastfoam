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
    Base class for a scalar equation with a single input

\*---------------------------------------------------------------------------*/

#ifndef equation_H
#define equation_H

#include "UnivariateEquation.H"
#include "ScalarEquation.H"
#include "univariateEquation.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class equation Declaration
\*---------------------------------------------------------------------------*/

class equation
:
    public Equation<scalar, scalar>,
    public scalarEquation,
    public univariateEquation<scalar>
{
public:

    // Constructors

        //- Construct given bounds
        equation(const scalar lowerLimit, const scalar upperLimit);


    //- Destructor
    virtual ~equation();


    // Member Functions

        //- Return the number of input variables
        inline label nVar() const;

        //- Return the number of equations
        inline label nEqns() const;

        using Equation<scalar, scalar>::nDerivatives;

        //- Return the lower limit
        inline tmp<scalarField> lowerLimits() const;

        //- Return the upper limit
        inline tmp<scalarField> upperLimits() const;

        //- Return the finite difference step size
        inline tmp<scalarField> dX() const;

        //- Return the finite difference step size
        inline void setDX(const scalarField& newDx) const;

        //- Limit the value of x to the range of the equation
        using scalarEquation::limit;
        inline void limit(scalar& x) const;

        //- Calulate the function at a point
        inline void f
        (
            const scalar& x,
            const label li,
            scalar& fx
        ) const;

        //- Calulate the function at a set of points
        inline void f
        (
            const scalarField& x,
            const label li,
            scalar& fx
        ) const;

        //- Calulate the gradient (dfdx)
        inline void gradient
        (
            const scalarField& x,
            const label li,
            scalar& fx,
            scalarField& grad
        ) const;

        //- Check if the solution is bracketed
        bool containsRoot(const scalar y0, const scalar y1) const;

        //- Check if the solution is bracketed (compute bounds)
        bool containsRoot(const label li) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
