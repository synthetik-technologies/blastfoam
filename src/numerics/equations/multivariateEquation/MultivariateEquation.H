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
    Foam::MultivariateEquation

Description
    Abstract base class a multivariate equation. Uses a scalar field as
    the input, but can output a list of any type.

\*---------------------------------------------------------------------------*/

#ifndef MultivariateEquation_H
#define MultivariateEquation_H

#include "multivariateEquation.H"
#include "Equation.H"
#include "scalarField.H"
#include "RectangularMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MultivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MultivariateEquation
:
    public Equation<scalarField, Field<Type>>,
    public multivariateEquation<Type>
{
protected:

    // Protected member functions

        //- Create a Jacobian matrix given a list of dx
        void calculateJacobian
        (
            const scalarField& x0,
            const label li,
            const Field<Type>& f0,
            RectangularMatrix<Type>& J
        ) const;

public:

    // Constructors

        //- Construct null
        MultivariateEquation(const label nEqns);

        //- Construct with limits
        MultivariateEquation
        (
            const label nEqns,
            const scalarField& lowerLimits,
            const scalarField& upperLimits
        );


    //- Destructor
    virtual ~MultivariateEquation();


    // Member Functions

        //- Return the number of input variables
        virtual label nVar() const;

        //- Return the number of equations
        virtual label nEqns() const;

        //- Return the lower limit
        virtual tmp<scalarField> lowerLimits() const;

        //- Return the upper limit
        virtual tmp<scalarField> upperLimits() const;

        //- Return the finite difference step size
        virtual tmp<scalarField> dX() const;

        //- Return the finite difference step size
        virtual void setDX(const scalarField& newDx) const;

        //- Limit the value of x to the range of the equation
        virtual void limit(scalarField& x) const;

        //- Calculate the Jacobian of the equation
        virtual void jacobian
        (
            const scalarField& x,
            const label li,
            Field<Type>& fx,
            RectangularMatrix<Type>& J
        ) const;
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
