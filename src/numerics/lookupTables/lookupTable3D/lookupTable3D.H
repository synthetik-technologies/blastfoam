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

    //- Include definition of modifying functions
    #include "tableFuncs.H"

    modFuncType modFunc_;
    modFuncType invModFunc_;

    modFuncType modXFunc_;
    modFuncType invModXFunc_;

    modFuncType modYFunc_;
    modFuncType invModYFunc_;

    modFuncType modZFunc_;
    modFuncType invModZFunc_;

    //- Pointer to function to lookup indexes in the x direction
    findIndexFunc findXIndex_;

    //- Pointer to function to lookup indexes in the y direction
    findIndexFunc findYIndex_;

    //- Pointer to function to lookup indexes in the y direction
    findIndexFunc findZIndex_;

    //- Interpolation type
    interp3DFuncType interpFunc_;

    //- Data
    Field<Field<Field<Type>>> data_;

    //- Modified x field values
    Field<scalar> xModValues_;

    //- Modified y field values
    Field<scalar> yModValues_;

    //- Modified z field values
    Field<scalar> zModValues_;

    //- Stored real x values
    Field<scalar> xValues_;

    //- Stored real y values
    Field<scalar> yValues_;

    //- Stored real z values
    Field<scalar> zValues_;

    //- Current indexes (used to access location in list after lookup)
    mutable label i_;
    mutable label j_;
    mutable label k_;

    //- Current linear interpolation weights
    mutable scalar fx_;
    mutable scalar fy_;
    mutable scalar fz_;

    //- Check if spacing in a list is uniform
    bool checkUniform(const List<scalar>& xyz) const;


public:

    // Constructors

        //- Construct empty
        lookupTable3D();

        //- Construct from dictionary
        lookupTable3D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& zName,
            const word& name
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
            const word& interpolationScheme,
            const bool isReal = true
        );



    //- Destructor
    virtual ~lookupTable3D();

    //- Access to data

        //- Access current index (x)
        label i() const
        {
            return i_;
        }

        //- Access current index (y)
        label j() const
        {
            return j_;
        }

        //- Access current index (z)
        label k() const
        {
            return k_;
        }

        //- Access current interpolation weight (x)
        scalar fx() const
        {
            return fx_;
        }

        //- Access current interpolation weight (y)
        scalar fy() const
        {
            return fy_;
        }

        //- Access current interpolation weight (z)
        scalar fz() const
        {
            return fz_;
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
            return invModXFunc_(f);
        }

        //- Modify by modYType
        scalar modY(const scalar& f) const
        {
            return modYFunc_(f);
        }

        //- Modify by invModYType
        scalar invModY(const scalar& f) const
        {
            return invModYFunc_(f);
        }

        //- Modify by modYType
        scalar modZ(const scalar& f) const
        {
            return modZFunc_(f);
        }

        //- Modify by invModYType
        scalar invModZ(const scalar& f) const
        {
            return invModZFunc_(f);
        }

        //- Const access to real x values
        const Field<scalar>& x() const
        {
            return xValues_;
        }

        //- Const access to real y values
        const Field<scalar>& y() const
        {
            return yValues_;
        }

        //- Const access to real z values
        const Field<scalar>& z() const
        {
            return zValues_;
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
        const Field<Field<Field<Type>>>& data() const
        {
            return data_;
        }

        //- Return real data values
        tmp<Field<Field<Field<Type>>>> realData() const;


    // Member Functions

        //- Update indexes
        void update(const scalar x, const scalar y, const scalar z) const;

        //- Update x index
        void updateX(const scalar x) const;

        //- Update y index
        void updateY(const scalar y) const;

        //- Update z index
        void updateZ(const scalar z) const;

        //- Update indexes
        void updateWeights
        (
            const scalar x,
            const scalar y,
            const scalar z
        ) const;

        //- Update x weight
        void updateXWeight(const scalar x) const;

        //- Update y weight
        void updateYWeight(const scalar y) const;

        //- Update z weight
        void updateZWeight(const scalar z) const;

        //- Lookup value
        Type lookup(const scalar x, const scalar y, const scalar z) const;

        //- Return value from a given list from lower index and weight
        Type getValue
        (
            const label ijk,
            const scalar f,
            const List<Type>& xyz
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
            const word& interpolationScheme,
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
            const word& name
        );

        //- Read a set from a given dictionary
        void readComponent
        (
            const dictionary& parentDict,
            const word& name,
            Field<scalar>& values,
            Field<scalar>& modValues,
            modFuncType& modFunc,
            modFuncType& invModFunc,
            findIndexFunc& findIndex
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "lookupTable3D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
