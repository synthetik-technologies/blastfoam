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
    Table used to lookup vales given a 1D table

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

    //- Include definition of modifying functions
    #include "tableFuncs.H"

    modFuncType modFunc_;
    modFuncType invModFunc_;

    modFuncType modXFunc_;
    modFuncType invModXFunc_;

    //- List of x values in real space
    scalarField xValues_;

    //- List of x values in the given space
    scalarField xModValues_;

    //- Data
    scalarField data_;

    //- Read the table
    void readTable();

    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline void findIndex(const scalar& x, label& I, scalar& f) const;

    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline label bound(const scalar& f) const;



public:

    // Constructors
    lookupTable1D
    (
        const fileName& fileName,
        const word& mod,
        const word& xMod
    );

    //- Destructor
    virtual ~lookupTable1D();


    // Member Functions

        //- Lookup value
        scalar lookup(const scalar& x) const;

        //- Lookup X given f and y
        scalar reverseLookup(const scalar& f) const;

        //- Return first derivative
        scalar dFdX(const scalar& x) const;

        //- Return second derivative w.r.t. x
        scalar d2FdX2(const scalar& x) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
