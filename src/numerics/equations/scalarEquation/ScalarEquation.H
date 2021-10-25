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

#ifndef ScalarEquation_H
#define ScalarEquation_H

#include "Equation.H"
#include "scalarEquation.H"
#include "scalarField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class scalarEquation Declaration
\*---------------------------------------------------------------------------*/

class ScalarEquation
:
    public Equation<scalarField, scalar>,
    public scalarEquation
{
protected:
// Protected data

    //- Calculate the gradient using a finite difference
    void calculateGradient
    (
        const scalarField& x0,
        const label li,
        const scalar& fx0,
        scalarField& grad
    ) const;


public:

    // Constructors

        //- Construct given bounds
        ScalarEquation
        (
            const label nVar,
            const scalarField& lowerLimit,
            const scalarField& upperLimit
        );


    //- Destructor
    virtual ~ScalarEquation();


    // Member Functions

        //- Return the number of input variables
        inline label nVar() const;

        //- Return the number of equations
        inline label nEqns() const;

        //- Return the lower limit
        inline tmp<scalarField> lowerLimits() const;

        //- Return the upper limit
        inline tmp<scalarField> upperLimits() const;

        //- Return the finite difference step size
        inline tmp<scalarField> dX() const;

        //- Return the finite difference step size
        inline void setDX(const scalarField& newDx) const;

        //- Limit the value of x to the range of the equation
        inline void limit(scalarField& x) const;

        //- Calulate the gradient
        using scalarEquation::gradient;
        virtual void gradient
        (
            const scalarField& x,
            const label li,
            scalar& fx,
            scalarField& grad
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ScalarEquationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
