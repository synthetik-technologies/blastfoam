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
    Table used to lookup values given a 1D table

SourceFiles
    lookupTable1D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable1D_H
#define lookupTable1D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "scalarList.H"
#include "IOstreams.H"
#include "Switch.H"
#include "IOList.H"
#include "fileName.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lookupTable1D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lookupTable1D
{
protected:
// Protected data

    //- Include definition of modifying functions
    #include "tableFuncs.H"

    modFuncType modFunc_;
    modFuncType invModFunc_;

    modFuncType modXFunc_;
    modFuncType invModXFunc_;

    //- Interpolation type
    interp1DFuncType interpFunc_;

    //- List of x values in real space
    scalarList xValues_;

    //- List of x values in the given space
    scalarList xModValues_;

    //- Data
    List<Type> data_;

    //- Current index (used to access location in list after lookup)
    mutable label index_;

    //- Current interpolation weight
    mutable scalar f_;

    //- Read the table
    void readTable
    (
        const fileName& file,
        const string& delim,
        const label xi,
        const label yi
    );


public:

    // Constructors

        //- Empty constructor
        lookupTable1D();

        //- Construct from file
        lookupTable1D
        (
            const dictionary& dict,
            const word& xName,
            const word& name
        );

        //- Construct from data given list of x and data
        //  x can be in real or xMod space
        lookupTable1D
        (
            const scalarList& x,
            const List<Type>& data,
            const word& xMod,
            const word& mod,
            const word& interpolationScheme,
            const bool isReal = true
        );

        //- Construct from interpolator from only x values, no data is stored
        lookupTable1D
        (
            const scalarList& x,
            const word& xMod,
            const word& interpolationScheme,
            const bool isReal = true
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
            Type mod(const Type& f) const
            {
                return modFunc_(f);
            }

            //- Modify by invModType
            Type invMod(const Type& f) const
            {
                return invModFunc_(f);
            }

            //- Modify by modXType
            scalar modX(const scalar f) const
            {
                return modXFunc_(f);
            }

            //- Modify by invModXType
            scalar invModX(const scalar f) const
            {
                return modXFunc_(f);
            }

            //- Const access to real x values
            const scalarList& x() const
            {
                return xValues_;
            }

            //- Const access to modified x values
            const scalarList& xMod() const
            {
                return xModValues_;
            }

            //- Const access to modified data values
            const List<Type>& data() const
            {
                return data_;
            }

            //- Return the real data
            tmp<List<Type>> realData() const;


    // Public functions

        //- Update index and weight
        void update(const scalar x) const;

        //- Lookup value
        Type lookup(const scalar x) const;

        //- Interpolate value using the current index and weight
        template<class fType>
        fType interpolate(const fType&, const fType&) const;

        //- Interpolate value using the current index and weight
        template<template<class> class ListType, class fType>
        fType interpolate(const scalar, const ListType<fType>&) const;

        //- Return first derivative
        Type dFdX(const scalar x) const;

        //- Return second derivative w.r.t. x
        Type d2FdX2(const scalar x) const;

        //- Set data
        void set
        (
            const scalarList& x,
            const List<Type>& data,
            const word& mod,
            const word& xMod,
            const word& interpolationScheme,
            const bool isReal
        );

        //- Set x
        void setX
        (
            const scalarList& x,
            const word& xMod,
            const word& interpolationScheme,
            const bool isReal = true
        );

        //- Set data
        void setData
        (
            const List<Type>& x,
            const word& mod,
            const bool isReal = true
        );

    //- Read from a given dictionary
    virtual void read
    (
        const dictionary& dict,
        const word& xName,
        const word& yName
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lookupTable1D.C"
    #include "lookupTable1DTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
