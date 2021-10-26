/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::HalleyUnivariateRootSolver

Description
    Find root using the Halley method

SourceFiles
    HalleyUnivariateRootSolver.C

\*---------------------------------------------------------------------------*/

#ifndef HalleyUnivariateRootSolver_H
#define HalleyUnivariateRootSolver_H

#include "univariateRootSolver.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HalleyUnivariateRootSolver Declaration
\*---------------------------------------------------------------------------*/

class HalleyUnivariateRootSolver
:
    public univariateRootSolver
{
public:

    //- Runtime type information
    TypeName("Halley");

    // Constructors

        //- Construct for given univariateRootSolver
        HalleyUnivariateRootSolver
        (
            const multivariateEquation<scalar>& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        HalleyUnivariateRootSolver(const univariateRootSolver&) = delete;


    //- Destructor
    virtual ~HalleyUnivariateRootSolver();


    // Member Functions

        //- Solve the root
        virtual scalar findRoot
        (
            const scalar x0,
            const scalar xLow,
            const scalar xHigh,
            const label li
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
