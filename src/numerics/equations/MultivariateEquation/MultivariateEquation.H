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
    Abstract base class a multivariate equation

\*---------------------------------------------------------------------------*/

#ifndef MultivariateEquation_H
#define MultivariateEquation_H

#include "scalarField.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class MultivariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MultivariateEquation
{
protected:

    // Protected data

        //- Lower limits
        scalarField lowerLimits_;

        //- Upper limits
        scalarField upperLimits_;


    // Protected member functions

        //- Make sure limits have been set
        void checkLimits() const;


public:

    // Constructors

        //- Construct null
        MultivariateEquation();

        //- Construct with limits
        MultivariateEquation
        (
            const scalarField& lowerLimits,
            const scalarField& upperLimits
        );


    //- Destructor
    virtual ~MultivariateEquation();


    // Member Functions

        //- Return the number of equations in the system
        virtual label nEqns() const = 0;

        //- Return the number of derivatives in the system
        virtual label nDerivatives() const = 0;

        //- Return the lower limits
        const scalarField& lower() const
        {
            return lowerLimits_;
        }

        //- Access the lower limits
        scalarField& lowerLimits()
        {
            return lowerLimits_;
        }

        //- Return the upper limits
        const scalarField upper() const
        {
            return upperLimits_;
        }

        //- Access the upper limits
        scalarField& upper()
        {
            return upperLimits_;
        }

        //- Check if the coordinates are with in the range
        bool checkBounds(const scalarField& xs) const;

        //- Return the function value
        tmp<Field<Type>> f
        (
            const scalarField& x,
            const label li
        ) const;

        //- Return the function value
        virtual void f
        (
            const scalarField& x,
            const label li,
            Field<Type>& fx
        ) const = 0;

        //- Calculate the Jacobian of the equation
        virtual void jacobian
        (
            const scalarField& x,
            const label li,
            scalarField& fx,
            SquareMatrix<Type>& dfdx
        ) const
        {
            NotImplemented;
        }
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
