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
    Foam::lookupTable2D

Description
    Table used to lookup values given a 2D table

SourceFiles
    lookupTable2D.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable2D_H
#define lookupTable2D_H

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
                           Class lookupTable2D Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class lookupTable2D
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

    //- Pointer to function to lookup indexes in the x direction
    findIndexFunc findXIndex_;

    //- Pointer to function to lookup indexes in the y direction
    findIndexFunc findYIndex_;

    //- Interpolation type
    interp2DFuncType interpFunc_;

    //- Data
    Field<Field<Type>> data_;

    //- Modified x field values
    Field<scalar> xModValues_;

    //- Modified y field values
    Field<scalar> yModValues_;

    //- Stored real x values
    Field<scalar> xValues_;

    //- Stored real y values
    Field<scalar> yValues_;

    //- Current indexes (used to access location in list after lookup)
    mutable label i_;
    mutable label j_;

    //- Current linear interpolation weights
    mutable scalar fx_;
    mutable scalar fy_;


    //- Read the table
    void readTable
    (
        const fileName& file,
        const string& delim,
        Field<Field<Type>>& data
    );

    //- Check if spacing in a list is uniform
    bool checkUniform(const List<scalar>& xy) const;


public:

    // Constructors

        //- Construct empty
        lookupTable2D();

        //- Construct from dictionary
        lookupTable2D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& name
        );

        //- Construct from data
        lookupTable2D
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<Field<Type>>& data,
            const word& xMod,
            const word& yMod,
            const word& mod,
            const word& interpolationScheme,
            const bool isReal = true
        );



    //- Destructor
    virtual ~lookupTable2D();

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

        //- Modify by modYType
        scalar modY(const scalar& f) const
        {
            return modYFunc_(f);
        }

        //- Modify by invModYType
        scalar invModY(const scalar& f) const
        {
            return modYFunc_(f);
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

        //- Const access to modified data values
        const Field<Field<Type>>& data() const
        {
            return data_;
        }

        //- Return real data values
        tmp<Field<Field<Type>>> realData() const;


    // Member Functions

        //- Update indexes
        void update(const scalar x, const scalar y) const;

        //- Update x index
        void updateX(const scalar x) const;

        //- Update y index
        void updateY(const scalar y) const;

        //- Update x weight
        void updateXWeight(const scalar x) const;

        //- Update y weight
        void updateYWeight(const scalar y) const;

        //- Lookup value
        Type lookup(const scalar x, const scalar y) const;

        //- Return first derivative w.r.t. x
        Type dFdX(const scalar x, const scalar y) const;

        //- Return first derivative w.r.t y
        Type dFdY(const scalar x, const scalar y) const;

        //- Return second derivative w.r.t. x
        Type d2FdX2(const scalar x, const scalar y) const;

        //- Return mixed second derivative
        Type d2FdXdY(const scalar x, const scalar y) const;

        //- Return second derivative w.r.t. y
        Type d2FdY2(const scalar x, const scalar y) const;

        //- Return value from a given list from lower index and weight
        Type getValue(const label ij, const scalar f, const List<Type>& xy) const;


    // Edit functions

        //- Set data, x, and y
        virtual void set
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<Field<Type>>& data,
            const bool isReal = true
        );
        virtual void set
        (
            const Field<scalar>& x,
            const Field<scalar>& y,
            const Field<Field<Type>>& data,
            const word& xMod,
            const word& yMod,
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

        //- Set data
        virtual void setData
        (
            const Field<Field<Type>>& data,
            const bool isReal = true
        );
        virtual void setData
        (
            const Field<Field<Type>>& data,
            const word& mod,
            const bool isReal = true
        );


        //- Read from a given dictionary
        virtual void read
        (
            const dictionary& dict,
            const word& xName,
            const word& yName,
            const word& name
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "lookupTable2D.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
