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
    Foam::scalarMultivariateRootSystem

Description
    Base class for a system of scalar equations

\*---------------------------------------------------------------------------*/

#ifndef scalarMultivariateEquation_H
#define scalarMultivariateEquation_H

#include "MultivariateEquation.H"
#include "scalarField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class scalarMultivariateEquation Declaration
\*---------------------------------------------------------------------------*/

class scalarMultivariateEquation
:
    public MultivariateEquation<scalar>
{
public:

    // Constructors

        //- Construct with limits
        scalarMultivariateEquation
        (
            const label nEqns,
            const scalarField& upperLimits,
            const scalarField& lowerLimits
        );


    //- Destructor
    virtual ~scalarMultivariateEquation();


    // Member Functions

        //- Check if the solution is bracketed
        bool containsRoot
        (
            const scalarField& y0s,
            const scalarField& y1s
        ) const;

        //- Check if the solution is bracketed (compute bounds)
        bool containsRoot(const label li) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
