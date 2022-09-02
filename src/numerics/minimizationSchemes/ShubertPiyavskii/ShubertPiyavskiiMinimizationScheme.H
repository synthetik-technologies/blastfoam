/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::ShubertPiyavskiiMinimizationScheme

Description
    Find root using the ShubertPiyavskii method

SourceFiles
    ShubertPiyavskiiMinimizationScheme.C

\*---------------------------------------------------------------------------*/

#ifndef ShubertPiyavskiiMinimizationScheme_H
#define ShubertPiyavskiiMinimizationScheme_H

#include "univariateMinimizationScheme.H"
#include "vector2D.H"
#include "typeInfo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ShubertPiyavskiiMinimizationScheme Declaration
\*---------------------------------------------------------------------------*/

class ShubertPiyavskiiMinimizationScheme
:
    public univariateMinimizationScheme
{
// Private data

    //- Lipschits constant
    scalar l_;

    //- Return intersection of two lines
    vector2D intersection(const vector2D&, const vector2D&) const;


public:

    //- Runtime type information
    TypeName("ShubertPiyavskii");

    // Constructors

        //- Construct for given minimizationScheme
        ShubertPiyavskiiMinimizationScheme
        (
            const scalarUnivariateEquation& eqn,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        ShubertPiyavskiiMinimizationScheme(const minimizationScheme&) = delete;


    //- Destructor
    virtual ~ShubertPiyavskiiMinimizationScheme()
    {}


    // Member Functions

        //- Minimize the function
        virtual scalar minimize
        (
            const scalar x0,
            const scalar x1,
            const scalar x2,
            const label li
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ShubertPiyavskiiMinimizationScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
