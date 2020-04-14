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
    Foam::equationOfStates::table

Description
    Table used to lookup equation of state properties

SourceFiles
    table.C

\*---------------------------------------------------------------------------*/

#ifndef lookupTable_H
#define lookupTable_H

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
                           Class lookupTable Declaration
\*---------------------------------------------------------------------------*/

class lookupTable
{
// Private data

    enum modType
    {
        none,
        log10,
        ln,
        exp
    };

    //- File name
    const fileName file_;

    //- Function for data set
    modType mod_;

    //- First direction function
    modType xMod_;

    //- Second direction modifier
    modType yMod_;

    label nx_;
    label ny_;

    scalar xMin_;
    scalar dx_;
    scalar yMin_;
    scalar dy_;

    //- Data
    Field<scalarField> data_;

    scalar getValue(const scalar& i, const scalar& minx, const scalar& dx) const
    {
        return minx + i*dx;
    }

    modType getModType(const word& type) const;

    //- Read the next value from the split string
    scalar readValue(const List<string>&) const;

    //- modify the variable based on the saved scheme
    scalar modVar(const modType& type, const scalar& xy) const;

    //- modify the variable based on the saved scheme
    scalar invModVar(const modType& type, const scalar& xy) const;

    //- Read the table
    void readTable();

    //- Find bottom of interpolation region, return index and weight between i and i+1
    inline void findIndex
    (
        const word& name,
        const scalar& xy,
        const scalar& xyMin,
        const scalar& dxy,
        const label nxy,
        const modType& mod,
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
        lookupTable
        (
            const word& fileName,
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
    virtual ~lookupTable();


    // Member Functions

        scalar lookup(const scalar& x, const scalar& y) const;

        scalar reverseLookupX(const scalar& f, const scalar& y) const;

        scalar reverseLookupY(const scalar& f, const scalar& x) const;

        scalar dFdX(const scalar& x, const scalar& y) const;

        scalar dFdY(const scalar& x, const scalar& y) const;

        scalar d2FdX2(const scalar& x, const scalar& y) const;

        scalar d2FdXdY(const scalar& x, const scalar& y) const;

        scalar d2FdY2(const scalar& x, const scalar& y) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
