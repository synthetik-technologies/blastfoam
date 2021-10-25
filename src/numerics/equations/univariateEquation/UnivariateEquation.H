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
    Foam::UnivariateEquation

Description
    Abstract base class a Univariate equation

\*---------------------------------------------------------------------------*/

#ifndef UnivariateEquation_H
#define UnivariateEquation_H

#include "Equation.H"
#include "univariateEquation.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UnivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UnivariateEquation
:
    public Equation<scalar, Type>,
    public univariateEquation<Type>
{
public:

    // Constructors

        //- Construct with limits
        UnivariateEquation
        (
            const label nVar,
            const scalar& lowerLimit,
            const scalar& upperLimit
        );


    //- Destructor
    virtual ~UnivariateEquation();


    // Member Functions

        //- Calculate the first derivative of the equation
        virtual void f
        (
            const scalar& x,
            const label li,
            Type& fx
        ) const
        {
            fx = this->fx(x, li);
        }
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
