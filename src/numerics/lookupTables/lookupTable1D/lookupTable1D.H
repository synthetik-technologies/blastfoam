/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
#include "Field.H"
#include "IOstreams.H"
#include "Switch.H"
#include "IOList.H"
#include "fileName.H"

#include "Modifier.H"
#include "indexing.H"
#include "interpolationWeights1D.H"

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

    //- Modifiers
    autoPtr<Modifier<Type>> mod_;
    autoPtr<Modifier<scalar>> modX_;

    //- Data
    Field<Type> data_;

    //- List of x values in the given space
    Field<scalar> xModValues_;

    //- Lookup indexes
    autoPtr<indexer> indexing_;

    //- Interpolation scheme
    autoPtr<interpolationWeight1D> interpolator_;

    //- Real data
    Field<Type>* realDataPtr_;

    //- List of x values in real space
    Field<scalar>* xValuesPtr_;

    //- Current index (used to access location in list after lookup)
    mutable label index_;

    //- List of interpolation indices
    mutable DynamicList<label> indices_;

    //- Current interpolation weight
    mutable DynamicList<scalar> weights_;


    // Protected member functions

        // Return the real data
        inline const Field<scalar>& fValues() const
        {
            return *realDataPtr_;
        }

        // Access the real data values
        inline Field<scalar>& fValues()
        {
            return *realDataPtr_;
        }

        // Return the real x values
        inline const Field<scalar>& xValues() const
        {
            return *xValuesPtr_;
        }

        // Access the real x values
        inline Field<scalar>& xValues()
        {
            return *xValuesPtr_;
        }


public:

    // Constructors

        //- Empty constructor
        lookupTable1D();

        //- Copy constructor
        lookupTable1D(const lookupTable1D&);

        //- Construct from file
        lookupTable1D
        (
            const dictionary& dict,
            const word& xName,
            const word& name,
            bool canRead = true
        );

        //- Construct from data given list of x and data
        //  x can be in real or xMod space
        lookupTable1D
        (
            const List<scalar>& x,
            const List<Type>& data,
            const word& xMod = "none",
            const word& mod = "none",
            const word& interpolationScheme  = "linearClamp",
            const bool isReal = true
        );

        //- Construct from interpolator from only x values, no data is stored
        lookupTable1D
        (
            const List<scalar>& x,
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

            //- Access current interpolation indices
            const labelList& indices() const
            {
                return indices_;
            }

            //- Access current interpolation indices
            const scalarList& weights() const
            {
                return weights_;
            }

            //- Return the modifier
            const Modifier<Type>& mod() const
            {
                return mod_();
            }

            //- Return the x-modifier
            const Modifier<scalar>& modX() const
            {
                return modX_();
            }

            //- Const access to real x values
            const Field<scalar>& x() const
            {
                return xValues();
            }

            //- Const access to modified x values
            const Field<scalar>& xMod() const
            {
                return xModValues_;
            }

            //- Const access to modified data values
            const Field<Type>& fMod() const
            {
                return data_;
            }

            //- Return the real data
            const Field<Type>& f() const
            {
                return fValues();
            }

            //- Return the interpolation scheme type name
            const interpolationWeight1D& interpolator() const
            {
                return interpolator_();
            }


    // Public functions

        //- Update index
        void updateIndex(const scalar x) const;

        //- Update index and weight
        void update(const scalar x) const;

        //- Lookup value
        Type lookup(const scalar x) const;

        //- Lookup x given f
        scalar reverseLookup(const Type& fin) const
        {
            NotImplemented;
            return 0.0;
        }

        //- Interpolate value using the current index and weight
        template<template<class> class ListType, class fType>
        fType interpolate(const ListType<fType>&) const;

        //- Interpolate value using the current index and weight
        template<template<class> class ListType, class fType>
        fType interpolate(const scalar, const ListType<fType>&) const;

        //- Return first derivative
        Type dFdX(const scalar x) const;

        //- Return second derivative w.r.t. x
        Type d2FdX2(const scalar x) const;

        //- Set x and data
        void set
        (
            const List<scalar>& x,
            const List<Type>& data,
            const bool isReal
        );
        void set
        (
            const List<scalar>& x,
            const List<Type>& data,
            const word& mod,
            const word& xMod,
            const word& interpolationScheme,
            const bool isReal
        );

        //- Set x
        void setX(const List<scalar>& x, const bool isReal = true);
        void setX
        (
            const List<scalar>& x,
            const word& xMod,
            const bool isReal = true
        );

        //- Set data
        void setData(const List<Type>& data, const bool isReal = true);
        void setData
        (
            const List<Type>& data,
            const word& mod,
            const bool isReal = true
        );

    //- Read from a given dictionary
    void read
    (
        const dictionary& dict,
        const word& xName,
        const word& yName,
        const bool canRead = true
    );
};


// Specilizations of scalar functions
template<>
scalar lookupTable1D<scalar>::reverseLookup(const scalar& fin) const;


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
