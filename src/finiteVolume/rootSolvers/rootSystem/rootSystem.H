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
    Foam::rootSystem

Description
    Abstract base class for the systems of ordinary differential equations.

\*---------------------------------------------------------------------------*/

#ifndef rootSystem_H
#define rootSystem_H

#include "scalarField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class rootSystem Declaration
\*---------------------------------------------------------------------------*/

class rootSystem
{
protected:

    // Protected data

        //- Lower limit
        scalar xMin_;

        //- Upper limit
        scalar xMax_;


public:

    // Constructors

        //- Construct null
        rootSystem();

        //- Construct given bounds
        rootSystem(const scalar xMin, const scalar xMax);


    //- Destructor
    virtual ~rootSystem();


    // Member Functions

        //- Return the lower limit
        scalar lower() const
        {
            return xMin_;
        }

        //- Access the lower limit
        scalar& lower()
        {
            return xMin_;
        }

        //- Return the upper limit
        scalar upper() const
        {
            return xMax_;
        }

        //- Access the upper limit
        scalar& upper()
        {
            return xMax_;
        }

        //- Return the function value
        virtual scalar f(const scalar x, const label li) const = 0;

        //- Calculate the first derivative of the equation
        virtual scalar dfdx(const scalar x, const label li) const = 0;

        //- Calculate the second derivative of the equation
        virtual scalar d2fdx2(const scalar x, const label li) const = 0;

        //- Check if the solution is bracketed
        void checkConditions(const scalar y0, const scalar y1) const;

        //- Check if the solution is bracketed (compute bounds)
        void checkConditions(const label li) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
