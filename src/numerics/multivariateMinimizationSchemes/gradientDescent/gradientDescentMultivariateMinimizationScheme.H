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
    Foam::gradientDescentMinimizationScheme

Description
    Find root using the gradientDescent method

SourceFiles
    gradientDescentMultivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef gradientDescentMultivariateMinimizationScheme_H
#define gradientDescentMultivariateMinimizationScheme_H

#include "multivariateMinimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class gradientDescentMultivariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class gradientDescentMultivariateMinimizationScheme
:
    public multivariateMinimizationScheme
{
public:

    //- Runtime type information
    TypeName("gradientDescent");

    // Constructors

        //- Construct for given minimizationScheme
        gradientDescentMultivariateMinimizationScheme
        (
            const scalarMultivariateEquation& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        gradientDescentMultivariateMinimizationScheme
        (
            const gradientDescentMultivariateMinimizationScheme&
        ) = delete;


    //- Destructor
    virtual ~gradientDescentMultivariateMinimizationScheme()
    {}


    // Member Functions

        //- Solve the root
        virtual tmp<scalarField> minimize
        (
            const scalarList& x0,
            const scalarList& x1,
            const scalarList& x2,
            const label li
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
