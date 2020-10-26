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

    //- Include definition of modifying functions
    #include "tableFuncs.H"

    word modType_;
    modFuncType modFunc_;
    modFuncType invModFunc_;

    word modXType_;
    modFuncType modXFunc_;
    modFuncType invModXFunc_;

    word modYType_;
    modFuncType modYFunc_;
    modFuncType invModYFunc_;

    //- Is the provided data uniform
    Switch uniformX_;
    Switch uniformY_;

    //- Number of entries for x
    label nx_;

    //- Number of entries for y
    label ny_;

    //- Data
    Field<scalarField> data_;

    //- Modified x field values
    scalarField xMod_;

    //- Modified y field values
    scalarField yMod_;

    //- Stored real x values
    scalarField x_;

    //- Stored real y values
    scalarField y_;

    //- Read the next value from the split string
    scalar readValue(const List<string>&) const;

    //- Read the table
    void readTable(const fileName& file, Field<scalarField>& data);

    //- Pointer to function to lookup indexes in the x direction
    void (*findXIndex_)(const scalar&, const scalarField&, label&, scalar&);

    //- Pointer to function to lookup indexes in the y direction
    void (*findYIndex_)(const scalar&, const scalarField&, label&, scalar&);

    //- Lookup based on uniform indexing
    static void findUniformIndexes
    (
        const scalar& xy,
        const scalarField& XY,
        label& IJ,
        scalar& f
    );

    //- Lookup based on non uniform indexing
    static void findNonuniformIndexes
    (
        const scalar& xy,
        const scalarField& XY,
        label& IJ,
        scalar& f
    );


    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline label boundi
    (
        const scalar& f,
        const Field<scalarField>& data,
        const label j
    ) const;

    //- Find bottom of interpolation region, return index and weight between j and j+1
    inline label boundj
    (
        const scalar& f,
        const Field<scalarField>& data,
        const label i
    ) const;

    //- Find bottom of interpolation region, return indexes and weights between i, i+1, j, and j+1
    inline void boundij
    (
        const scalar& f,
        const Field<scalarField>& data,
        label& i,
        label& j
    ) const;

    //- Return value from a given list from lower index and weight
    scalar getValue(const label ij, const scalar& f, const scalarField& xy) const;

    //- Check if spacing in a list is uniform
    bool checkUniform(const scalarField& xy) const;

public:

    // Constructors

        //- Construct from dictionary
        lookupTable2D
        (
            const dictionary& dict,
            const word& xName,
            const word& yName
        );

        //- Construct from file
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

        //- Construct from data
        lookupTable2D
        (
            const Field<scalarField>& data,
            const scalarField& x,
            const scalarField& y,
            const word& mod,
            const word& xMod,
            const word& yMod,
            const bool modified
        );



    //- Destructor
    virtual ~lookupTable2D();

    //- Access to data

        //- Return name of modifier
        const word& modType() const
        {
            return modType_;
        }

        //- Return name of x modifier
        const word& modXType() const
        {
            return modXType_;
        }

        //- Return name of y modifier
        const word& modYType() const
        {
            return modYType_;
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
        const scalarField& x() const
        {
            return x_;
        }

        //- Non-const access to real x values
        scalarField& x()
        {
            return x_;
        }

        //- Const access to real y values
        const scalarField& y() const
        {
            return y_;
        }

        //- Non-const access to real y values
        scalarField& y()
        {
            return y_;
        }

        //- Const access to modified x values
        const scalarField& xMod() const
        {
            return xMod_;
        }

        //- Non-const access to modified x values
        scalarField& xMod()
        {
            return xMod_;
        }

        //- Const access to modified y values
        const scalarField& yMod() const
        {
            return yMod_;
        }

        //- Non-const access to modified y values
        scalarField& yMod()
        {
            return yMod_;
        }

        //- Const access to modified data values
        const Field<scalarField>& data() const
        {
            return data_;
        }

        //- Non-const access to modified data values
        Field<scalarField>& data()
        {
            return data_;
        }


    // Member Functions

        //- Return real data values
        tmp<Field<scalarField>> realData() const;

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
