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
    Foam::lookupTable3D

Description
    Table used to lookup values given a 3D table

SourceFiles
    lookupTable3D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable3D_H
#define lookupTable3D_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "Field.H"
#include "IOstreams.H"
#include "Switch.H"
#include "Field.H"
#include "fileName.H"
#include "labelVector.H"

#include "Modifier.H"
#include "interpolationWeights1D.H"
#include "indexing.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class lookupTable3D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lookupTable3D
{
protected:
// Protected data

    //- Modifiers
    autoPtr<Modifier<Type>> mod_;
    autoPtr<Modifier<scalar>> modX_;
    autoPtr<Modifier<scalar>> modY_;
    autoPtr<Modifier<scalar>> modZ_;

    //- Data
    Field<Field<Field<Type>>> data_;

    //- Modified x field values
    Field<scalar> xModValues_;

    //- Modified y field values
    Field<scalar> yModValues_;

    //- Modified z field values
    Field<scalar> zModValues_;

    //- Lookup indexes in the x direction
    autoPtr<indexer> xIndexing_;

    //- Lookup indexes in the y direction
    autoPtr<indexer> yIndexing_;

    //- Lookup indexes in the z direction
    autoPtr<indexer> zIndexing_;

    //- X-interpolator
    autoPtr<interpolationWeight1D> xInterpolator_;

    //- Y-interpolator
    autoPtr<interpolationWeight1D> yInterpolator_;

    //- Z-interpolator
    autoPtr<interpolationWeight1D> zInterpolator_;

    //- Stored real data
    Field<Field<Field<Type>>>* realDataPtr_;

    //- Stored real x values
    Field<scalar>* xValuesPtr_;

    //- Stored real y values
    Field<scalar>* yValuesPtr_;

    //- Stored real z values
    Field<scalar>* zValuesPtr_;

    //- Current indexes (used to access location in list after lookup)
    mutable labelVector ijk_;
    mutable List<labelVector> indices_;

    //- Current interpolation weights
    mutable scalarList weights_;

    //- Check if spacing in a list is uniform
    bool checkUniform(const List<scalar>& xyz) const;


    // Protected member functions

        //- Return the real data
        inline const Field<Field<Field<scalar>>>& fValues() const
        {
            return *realDataPtr_;
        }

        //- Access the real data values
        inline Field<Field<Field<scalar>>>& fValues()
        {
            return *realDataPtr_;
        }

        //- Return the real x values
        inline const Field<scalar>& xValues() const
        {
            return *xValuesPtr_;
        }

        //- Access the real x values
        inline Field<scalar>& xValues()
        {
            return *xValuesPtr_;
        }

        //- Return the real y values
        inline const Field<scalar>& yValues() const
        {
            return *yValuesPtr_;
        }

        //- Access the real y values
        inline Field<scalar>& yValues()
        {
            return *yValuesPtr_;
        }

        //- Return the real z values
        inline const Field<scalar>& zValues() const
        {
            return *zValuesPtr_;
        }

        //- Access the real z values
        inline Field<scalar>& zValues()
        {
            return *zValuesPtr_;
        }


public:

    // Constructors

        //- Construct empty
        lookupTable3D();

        //- Copy constructor
        lookupTable3D(const lookupTable3D&);

        //- Construct from dictionary
        lookupTable3D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& zName,
            const word& name,
            const bool canRead = true
        );

        //- Construct from data
        lookupTable3D
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<scalar>& z,
            const Field<Field<Field<Type>>>& data,
            const word& xMod,
            const word& yMod,
            const word& zMod,
            const word& mod,
            const word& xInterpolationScheme,
            const word& yInterpolationScheme,
            const word& zInterpolationScheme,
            const bool isReal = true
        );



    //- Destructor
    virtual ~lookupTable3D();

    //- Access to data

        //- Access current index (x)
        label i() const
        {
            return ijk_.x();
        }

        //- Access current index (y)
        label j() const
        {
            return ijk_.y();
        }

        //- Access current index (z)
        label k() const
        {
            return ijk_.z();
        }

        //- Access current interpolation indices
        const List<labelVector>& indices() const
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

        //- Return the y-modifier
        const Modifier<scalar>& modY() const
        {
            return modY_();
        }

        //- Return the z-modifier
        const Modifier<scalar>& modZ() const
        {
            return modZ_();
        }

        //- Const access to real x values
        const Field<scalar>& x() const
        {
            return xValues();
        }

        //- Const access to real y values
        const Field<scalar>& y() const
        {
            return yValues();
        }

        //- Const access to real z values
        const Field<scalar>& z() const
        {
            return zValues();
        }

        //- Const access to modified x values
        const Field<scalar>& xMod() const
        {
            return xModValues_;
        }

        //- Const access to modified y values
        const Field<scalar>& yMod() const
        {
            return yModValues_;
        }

        //- Const access to modified z values
        const Field<scalar>& zMod() const
        {
            return zModValues_;
        }

        //- Const access to modified data values
        const Field<Field<Field<Type>>>& fMod() const
        {
            return data_;
        }

        //- Return real data values
        const Field<Field<Field<Type>>>& f() const
        {
            return fValues();
        }

        //- Return the interpolation scheme for x
        const interpolationWeight1D& xInterpolator() const
        {
            return xInterpolator_();
        }

        //- Return the interpolation scheme for y
        const interpolationWeight1D& yInterpolator() const
        {
            return yInterpolator_();
        }

        //- Return the interpolation scheme for z
        const interpolationWeight1D& zInterpolator() const
        {
            return zInterpolator_();
        }


    // Member Functions

        //- Update indexes
        void updateIndex(const scalar x, const scalar y, const scalar z) const;

        //- Update indexes and weights
        void update(const scalar x, const scalar y, const scalar z) const;

        //- Lookup value
        Type lookup(const scalar x, const scalar y, const scalar z) const;

        //- Return value from a given list from lower index and weight
        Type getValue
        (
            const label ijk,
            const scalar f,
            const List<Type>& xyz
        ) const;

        //- Interpolate value using the current index and weight
        template
        <
            template<class> class ListType1,
            template<class> class ListType2,
            template<class> class ListType3,
            class fType
        >
        fType interpolate
        (
            const ListType1<ListType2<ListType3<fType>>>&
        ) const;

        //- Interpolate value using the current index and weight
        template
        <
            template<class> class ListType1,
            template<class> class ListType2,
            template<class> class ListType3,
            class fType
        >
        fType interpolate
        (
            const scalar,
            const scalar,
            const scalar,
            const ListType1<ListType2<ListType3<fType>>>&
        ) const;


    // Edit functions

        //- Set data, x, and y
        virtual void set
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<scalar>& z,
            const Field<Field<Field<Type>>>& data,
            const bool isReal = true
        );
        virtual void set
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<scalar>& z,
            const Field<Field<Field<Type>>>& data,
            const word& xMod,
            const word& yMod,
            const word& zMod,
            const word& mod,
            const word& xInterpolationScheme,
            const word& yInterpolationScheme,
            const word& zInterpolationScheme,
            const bool isReal = true
        );

        //- Set x
        virtual void setX
        (
            const Field<scalar>& x,
            const bool isReal = true
        );
        virtual void setX
        (
            const Field<scalar>& x,
            const word& xMod,
            const bool isReal = true
        );

        //- Set y
        virtual void setY
        (
            const Field<scalar>& y,
            const bool isReal = true
        );
        virtual void setY
        (
            const Field<scalar>& y,
            const word& yMod,
            const bool isReal = true
        );

        //- Set z
        virtual void setZ
        (
            const Field<scalar>& z,
            const bool isReal = true
        );
        virtual void setZ
        (
            const Field<scalar>& z,
            const word& zMod,
            const bool isReal = true
        );

        //- Set data
        virtual void setData
        (
            const Field<Field<Field<Type>>>& data,
            const bool isReal = true
        );
        virtual void setData
        (
            const Field<Field<Field<Type>>>& data,
            const word& mod,
            const bool isReal = true
        );


        //- Read from a given dictionary
        virtual void read
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& zName,
            const word& name,
            const bool canRead = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "lookupTable3D.C"
    #include "lookupTable3DTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
