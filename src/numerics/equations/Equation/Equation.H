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
    Foam::Equation

Description
    Abstract base class for the an equation.

\*---------------------------------------------------------------------------*/

#ifndef Equation_H
#define Equation_H

#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Equation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Equation
{
protected:

    // Protected data

        //- Lower limit
        scalar lowerLimit_;

        //- Upper limit
        scalar upperLimit_;


public:

    // Constructors

        //- Construct null
        Equation();

        //- Construct given bounds
        Equation(const Type lowerLimit, const Type upperLimit);


    //- Destructor
    virtual ~Equation();


    // Member Functions

        //- Return the lower limit
        scalar lower() const
        {
            return lowerLimit_;
        }

        //- Access the lower limit
        scalar& lower()
        {
            return lowerLimit_;
        }

        //- Return the upper limit
        scalar upper() const
        {
            return upperLimit_;
        }

        //- Access the upper limit
        scalar& upper()
        {
            return upperLimit_;
        }

        //- Check if derivatives have been implemented
        static bool checkFirstDerivative(const Equation<Type>& eqn);
        static bool checkSecondDerivative(const Equation<Type>& eqn);
        static bool checkThirdDerivative(const Equation<Type>& eqn);
        static bool checkFourthDerivative(const Equation<Type>& eqn);

        //- Return the number of derivatives
        label nDerivatives() const;

        //- Check if the coordinates are with in the range
        bool checkBounds(const scalar xs) const;

        //- Return the function value
        virtual Type f(const scalar x, const label li) const = 0;

        //- Calculate the first derivative of the equation
        virtual Type dfdx(const scalar x, const label li) const
        {
            NotImplemented
            return pTraits<Type>::zero;
        }

        //- Calculate the second derivative of the equation
        virtual Type d2fdx2(const scalar x, const label li) const
        {
            NotImplemented
            return pTraits<Type>::zero;
        }

        //- Calculate the third derivative of the equation
        virtual Type d3fdx3(const scalar x, const label li) const
        {
            NotImplemented;
            return pTraits<Type>::zero;
        }

        //- Calculate the fourth derivative of the equation
        virtual Type d4fdx4(const scalar x, const label li) const
        {
            NotImplemented;
            return pTraits<Type>::zero;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "Equation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
