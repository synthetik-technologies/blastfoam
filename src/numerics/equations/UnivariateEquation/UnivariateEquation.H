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
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UnivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class EquationType>
class UnivariateEquation
:
    public EquationType
{
public:

    typedef typename EquationType::inType inType;
    typedef typename EquationType::outType outType;

    // Constructors
    
        //- Construct with limits
        UnivariateEquation
        (
            const label nVar,
            const inType& lowerLimit,
            const inType& upperLimit
        );


    //- Destructor
    virtual ~UnivariateEquation();


    // Member Functions

        //- Calculate the first derivative of the equation
        using EquationType::equationType::f;
        virtual void f
        (
            const inType& x,
            const label li,
            outType& fx
        ) const
        {
            fx = this->fx(x, li);
        }

        virtual outType fx(const inType& x, const label li) const = 0;

        //- Calculate the first derivative of the equation
        virtual outType dfdx(const inType& x, const label li) const
        {
            NotImplemented
            return this->fx(x, li);
        }

        //- Calculate the second derivative of the equation
        virtual outType d2fdx2(const inType& x, const label li) const
        {
            NotImplemented
            return this->fx(x, li);
        }

        //- Calculate the third derivative of the equation
        virtual outType d3fdx3(const inType& x, const label li) const
        {
            NotImplemented;
            return this->fx(x, li);
        }

        //- Calculate the fourth derivative of the equation
        virtual outType d4fdx4(const inType& x, const label li) const
        {
            NotImplemented;
            return this->fx(x, li);
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
