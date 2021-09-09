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
    Foam::scalarLookupTable2D

Description
    Table used to lookup scalar values given a 2D table

SourceFiles
    scalarLookupTable2D.C

\*---------------------------------------------------------------------------*/

#ifndef scalarLookupTable2D_H
#define scalarLookupTable2D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lookupTable2D.H"
#include "ListOps.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scalarLookupTable2D Declaration
\*---------------------------------------------------------------------------*/

class scalarLookupTable2D
:
    public lookupTable2D<scalar>
{
// Private data

    //- Include definition of modifying functions
    #include "scalarTableFuncs.H"

    //- Find bottom of interpolation region, return index and weight between i and i+1
    labelList boundi(const scalar f) const;

    //- Find bottom of interpolation region, return index and weight between j and j+1
    labelList boundj(const scalar f) const;


public:

    // Constructors

        //- Construct empty
        scalarLookupTable2D();

        //- Construct from dictionary
        scalarLookupTable2D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& name
        );


        //- Construct from data
        scalarLookupTable2D
        (
            const scalarList& x,
            const scalarList& y,
            const List<List<scalar>>& data,
            const word& xMod,
            const word& yMod,
            const word& mod,
            const word& interpolationScheme,
            const bool isReal
        );



    //- Destructor
    virtual ~scalarLookupTable2D();

    // Member Functions

        //- Lookup X given f and y
        scalar reverseLookupX(const scalar f, const scalar y) const;

        //- Lookup y given f and x
        scalar reverseLookupY(const scalar f, const scalar x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
