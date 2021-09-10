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
    Foam::multivariateRootSystem

Description
    Abstract base class for the systems of ordinary differential equations.

\*---------------------------------------------------------------------------*/

#ifndef multivariateRootSystem_H
#define multivariateRootSystem_H

#include "scalarField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class multivariateRootSystem Declaration
\*---------------------------------------------------------------------------*/

class multivariateRootSystem
{
protected:

    // Protected data

        //- Lower limits
        scalarField xMins_;

        //- Upper limits
        scalarField xMaxs_;


    // Protected member functions

        //- Make sure limits have been set
        void checkLimits() const;

public:

    // Constructors

        //- Construct null
        multivariateRootSystem();

        //- Construct with limits
        multivariateRootSystem
        (
            const scalarField& xMins,
            const scalarField& xMaxs
        );


    //- Destructor
    virtual ~multivariateRootSystem();


    // Member Functions

        //- Return the number of equations in the system
        virtual label nEqns() const = 0;

        //- Return the lower limits
        const scalarField& lower() const
        {
            return xMins_;
        }

        //- Access the lower limits
        scalarField& lower()
        {
            return xMins_;
        }

        //- Return the upper limits
        const scalarField upper() const
        {
            return xMaxs_;
        }

        //- Access the upper limits
        scalarField& upper()
        {
            return xMaxs_;
        }

        //- Return the function value
        virtual void f
        (
            const scalarField& x,
            const label li,
            scalarField& fx
        ) const = 0;

        //- Calculate the first derivative of the equation
        virtual void jacobian
        (
            const scalarField& x,
            const label li,
            scalarField& fx,
            scalarSquareMatrix& dfdx
        ) const = 0;

        //- Check if the solution is bracketed
        void checkConditions
        (
            const scalarField& y0s,
            const scalarField& y1s
        ) const;

        //- Check if the solution is bracketed (compute bounds)
        void checkConditions(const label li) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
