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
    Foam::univariateEquation

Description
    Abstract base class an equation with a single scalar input, and any type
    as an output.

\*---------------------------------------------------------------------------*/

#ifndef univariateEquation_H
#define univariateEquation_H

#include "Equation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class univariateEquation Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class univariateEquation
{
public:

    // Constructors

        //- Construct with limits
        univariateEquation()
        {}


    //- Destructor
    virtual ~univariateEquation()
    {}


    // Member Functions

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


        //- Optional printing of stringified version of the equations
        virtual string printfx() const
        {
            NotImplemented;
            return string::null;
        }
        virtual string printdfdx() const
        {
            NotImplemented;
            return string::null;
        }
        virtual string printd2fdx2() const
        {
            NotImplemented;
            return string::null;
        }
        virtual string printd3fdx3() const
        {
            NotImplemented;
            return string::null;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
