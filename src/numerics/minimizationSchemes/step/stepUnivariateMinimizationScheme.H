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
    Foam::stepUnivariateMinimizationScheme

Description
    Fine the minimum of a function by stepping by a given increment. If the
    value at the next location is large, the step size is reduced by the
    given factor

SourceFiles
    stepUnivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef stepUnivariateMinimizationScheme_H
#define stepUnivariateMinimizationScheme_H

#include "univariateMinimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class stepUnivariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class stepUnivariateMinimizationScheme
:
    public univariateMinimizationScheme
{
// Private data

    //- Initial step size
    scalar dx_;

    //- Step size reduction
    scalar f_;


public:

    //- Runtime type information
    TypeName("step");

    // Constructors

        //- Construct for given univariateMinimizationScheme
        stepUnivariateMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

    //- Destructor
    virtual ~stepUnivariateMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        virtual scalar minimize
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const stepUnivariateMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
