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
    Foam::lookupTable1D

Description
    Table used to a one dimensional set of data

SourceFiles
    lookupTable1D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable1D_H
#define lookupTable1D_H

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
                           Class lookupTable1D Declaration
\*---------------------------------------------------------------------------*/

class lookupTable1D
{
// Private data

    //- File name
    const fileName file_;

    scalar (*modFunc_)(const scalar&);
    scalar (*invModFunc_)(const scalar&);

    scalar (*modXFunc_)(const scalar&);
    scalar (*invModXFunc_)(const scalar&);

    //- Data
    scalarField data_;

    //- Stored table lists in real space
    scalarField x_;
    scalarField realX_;

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
        const scalar& x,
        label& I,
        scalar& f
    ) const;

    inline static scalar noneS(const scalar& x)
    {
        return x;
    }

    inline static scalar lnS(const scalar& x)
    {
        return log(max(x, small));
    }

    inline static scalar expS(const scalar& x)
    {
        return exp(x);
    }

    inline static scalar log10S(const scalar& x)
    {
        return Foam::log10(max(x, small));
    }

    inline static scalar pow10S(const scalar& x)
    {
        return Foam::pow(10, x);
    }


public:

    // Constructors

        //- Multiphase constructor
        lookupTable1D
        (
            const fileName& fileName,
            const word& mod,
            const word& xMod
        );


    //- Destructor
    virtual ~lookupTable1D();


    // Member Functions

        //- Look up value
        scalar lookup(const scalar& x) const;

        //- Return f given x
        scalar reverseLookup(const scalar& f) const;

        //- Return first derivative
        scalar dFdX(const scalar& x) const;

        //- Return Second derivative
        scalar d2FdX2(const scalar& x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
