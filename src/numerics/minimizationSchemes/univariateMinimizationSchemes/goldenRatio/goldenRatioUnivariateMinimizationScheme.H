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
    Foam::goldenRatioUnivariateMinimizationScheme

Description
    Find root using the goldenRatio method

SourceFiles
    goldenRatioUnivariateMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef goldenRatioUnivariateMinimizationScheme_H
#define goldenRatioUnivariateMinimizationScheme_H

#include "univariateMinimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class goldenRatioUnivariateMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class goldenRatioUnivariateMinimizationScheme
:
    public univariateMinimizationScheme
{
    const static scalar goldenRatio;
    const static scalar invPhi;
    const static scalar invPhi2;

public:

    //- Runtime type information
    TypeName("goldenRatio");

    // Constructors

        //- Construct for given univariateUnivariateMinimizationScheme
        goldenRatioUnivariateMinimizationScheme
        (
            const scalarEquation& eqn,
            const dictionary& dict
        );


    //- Destructor
    virtual ~goldenRatioUnivariateMinimizationScheme()
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
        void operator=(const goldenRatioUnivariateMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
