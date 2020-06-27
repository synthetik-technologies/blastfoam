/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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
    Foam::lookupTable2D

Description
    Table used to lookup vales given a 2D table

SourceFiles
    lookupTable2D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable2D_H
#define lookupTable2D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "Field.H"
#include "fieldTypes.H"
#include "FieldFields.H"
#include "IOstreams.H"
#include "Switch.H"
#include "IOField.H"
#include "fileName.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lookupTable2D Declaration
\*---------------------------------------------------------------------------*/

class lookupTable2D
{
// Private data

    //- File name
    const fileName file_;

    //- Include definition of modifying functions
    #include "tableFuncs.H"

    modFuncType modFunc_;
    modFuncType invModFunc_;

    modFuncType modXFunc_;
    modFuncType invModXFunc_;

    modFuncType modYFunc_;
    modFuncType invModYFunc_;

    //- Number of entries for x
    label nx_;

    //- Number of entries for y
    label ny_;

    //- Minimum x value (in given space)
    scalar xMin_;

    //- X spacing (in given space)
    scalar dx_;

    //- Minimum y value (in given space)
    scalar yMin_;

    //- Y spacing (in given space)
    scalar dy_;

    //- Data
    Field<scalarField> data_;

    //- Stored table lists in real space
    scalarField x_;
    scalarField y_;

    //- Return x or y given the index
    inline scalar getValue(const scalar& i, const scalar& minx, const scalar& dx) const
    {
        return minx + i*dx;
    }

    //- Read the next value from the split string
    scalar readValue(const List<string>&) const;

    //- Read the table
    void readTable();

    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline void findIndex
    (
        const scalar& xy,
        const scalar& xyMin,
        const scalar& dxy,
        const label nxy,
        label& IJ,
        scalar& f
    ) const;



    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline label bound
    (
        const scalar& f,
        const label i,
        const scalar& yMin,
        const scalar& dy,
        const label ny,
        const bool ij
    ) const;


public:

    // Constructors

        //- Multiphase constructor
        lookupTable2D
        (
            const fileName& fileName,
            const word& mod,
            const word& xMod,
            const word& yMod,
            const label nx,
            const label ny,
            const scalar& xMin,
            const scalar& dx,
            const scalar& yMin,
            const scalar& dy
        );



    //- Destructor
    virtual ~lookupTable2D();


    // Member Functions

        //- Lookup value
        scalar lookup(const scalar& x, const scalar& y) const;

        //- Lookup X given f and y
        scalar reverseLookupX(const scalar& f, const scalar& y) const;

        //- Lookup y given f and x
        scalar reverseLookupY(const scalar& f, const scalar& x) const;

        //- Return first derivative w.r.t. x
        scalar dFdX(const scalar& x, const scalar& y) const;

        //- Return first derivative w.r.t y
        scalar dFdY(const scalar& x, const scalar& y) const;

        //- Return second derivative w.r.t. x
        scalar d2FdX2(const scalar& x, const scalar& y) const;

        //- Return mixed second derivative
        scalar d2FdXdY(const scalar& x, const scalar& y) const;

        //- Return second derivative w.r.t. y
        scalar d2FdY2(const scalar& x, const scalar& y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
