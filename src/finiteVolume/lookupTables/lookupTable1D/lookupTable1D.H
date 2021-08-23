/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

    //- Include definition of modifying functions
    #include "tableFuncs.H"

    modFuncType modFunc_;
    modFuncType invModFunc_;

    modFuncType modXFunc_;
    modFuncType invModXFunc_;

    //- Delimiter
    string delim_ = ",";

    //- Interpolation type
    interpFuncType interpFunc_;

    //- List of x values in real space
    scalarField xValues_;

    //- List of x values in the given space
    scalarField xModValues_;

    //- Data
    scalarField data_;

    //- Current index (used to access location in list after lookup)
    mutable label index_;

    //- Current interpolation weight
    mutable scalar f_;

    //- Read the table
    void readTable
    (
        const fileName& file,
        const label xi = 0,
        const label yi = 1,
        bool isReal = false
    );

    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline label bound(const scalar& f) const;



public:

    // Constructors

        //- Empty constructor
        lookupTable1D();

        //- Construct from file
        lookupTable1D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName
        );

        //- Construct from file
        lookupTable1D
        (
            const fileName& file,
            const word& mod,
            const word& xMod,
            const word& interpolationScheme = "linearClamp"
        );

        //- Construct from data given list of x and data
        //  x can be in real or xMod space
        lookupTable1D
        (
            const scalarField& x,
            const scalarField& data,
            const word& mod,
            const word& xMod,
            const word& interpolationScheme,
            const bool correct = false
        );

        //- Construct from interpolator from only x values, no data is stored
        lookupTable1D
        (
            const scalarField& x,
            const word& xMod,
            const word& interpolationScheme,
            const bool correct = false
        );


    //- Destructor
    virtual ~lookupTable1D();


    // Member Functions

        // Access data

            //- Access current index
            label index() const
            {
                return index_;
            }

            //- Access current interpolation weight
            scalar f() const
            {
                return f_;
            }

            //- Modify by modType
            scalar mod(const scalar& f) const
            {
                return modFunc_(f);
            }

            //- Modify by invModType
            scalar invMod(const scalar& f) const
            {
                return invModFunc_(f);
            }

            //- Modify by modXType
            scalar modX(const scalar& f) const
            {
                return modXFunc_(f);
            }

            //- Modify by invModXType
            scalar invModX(const scalar& f) const
            {
                return modXFunc_(f);
            }

            //- Const access to real x values
            const scalarField& x() const
            {
                return xValues_;
            }

            //- Const access to modified x values
            const scalarField& xMod() const
            {
                return xModValues_;
            }

            //- Non-const access to modified x values
            scalarField& xMod()
            {
                return xModValues_;
            }

            //- Const access to modified data values
            const scalarField& data() const
            {
                return data_;
            }

            //- Non-const access to modified data values
            scalarField& data()
            {
                return data_;
            }


    // Public functions

        //- Update index and weight
        void update(const scalar& x) const;

        //- Lookup value
        scalar lookup(const scalar& x) const;

        //- Interpolate value using the current index and weight
        template<class Type>
        Type interpolate(const Type&, const Type&) const;

        //- Interpolate value using the current index and weight
        template<template<class> class ListType, class Type>
        Type interpolate(const scalar, const ListType<Type>&) const;


        //- Lookup X given f and y
        scalar reverseLookup(const scalar& f) const;

        //- Return first derivative
        scalar dFdX(const scalar& x) const;

        //- Return second derivative w.r.t. x
        scalar d2FdX2(const scalar& x) const;

        //- Set data
        void set
        (
            const scalarField& x,
            const scalarField& data,
            const word& mod,
            const word& xMod,
            const word& interpolationScheme,
            const bool inReal
        );

        //- Set data
        void setX
        (
            const scalarField& x,
            const bool inReal = true
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lookupTable1DTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
