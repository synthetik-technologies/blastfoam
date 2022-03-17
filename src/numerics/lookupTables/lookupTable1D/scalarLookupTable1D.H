/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::scalarLookupTable1D

Description
    Table used to lookup scalar values given a 1D table

SourceFiles
    scalarLookupTable1D.C

\*---------------------------------------------------------------------------*/

#ifndef scalarLookupTable1D_H
#define scalarLookupTable1D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lookupTable1D.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class scalarLookupTable1D Declaration
\*---------------------------------------------------------------------------*/

class scalarLookupTable1D
:
    public lookupTable1D<scalar>
{
// Private data

public:

    // Constructors

        //- Empty constructor
        scalarLookupTable1D();

        //- Construct from file
        scalarLookupTable1D
        (
            const dictionary& dict,
            const word& xName,
            const word& name,
            const bool canRead = true
        );

        //- Construct from data given list of x and data
        //  x can be in real or xMod space
        scalarLookupTable1D
        (
            const Field<scalar>& x,
            const Field<scalar>& data,
            const word& mod,
            const word& xMod,
            const word& interpolationScheme,
            const bool correct = false
        );

        //- Construct from interpolator from only x values, no data is stored
        scalarLookupTable1D
        (
            const Field<scalar>& x,
            const word& xMod,
            const word& interpolationScheme,
            const bool correct = false
        );


    //- Destructor
    virtual ~scalarLookupTable1D();


    // Public functions

        //- Lookup X given f
        scalar reverseLookup(const scalar f) const;

        //- Read from a given dictionary
        virtual void read
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const bool canRead = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
