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
    Foam::goldenRatioMinimizationScheme

Description
    Find root using the goldenRatio method

SourceFiles
    goldenRatioMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef goldenRatioMinimizationScheme_H
#define goldenRatioMinimizationScheme_H

#include "minimizationScheme.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class goldenRatioMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class goldenRatioMinimizationScheme
:
    public minimizationScheme
{
    const static scalar goldenRatio;
    const static scalar invPhi;
    const static scalar invPhi2;

public:

    //- Runtime type information
    TypeName("goldenRatio");

    // Constructors

        //- Construct for given minimizationScheme
        goldenRatioMinimizationScheme(const scalarEquation& eqn, const dictionary& dict);

        //- Disallow default bitwise copy construction
        goldenRatioMinimizationScheme(const minimizationScheme&) = delete;


    //- Destructor
    virtual ~goldenRatioMinimizationScheme()
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
        void operator=(const goldenRatioMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
