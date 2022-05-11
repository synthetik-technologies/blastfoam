/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
    Foam::cellStencil

Description
    Basic class for a cell stencil

SourceFiles
    cellStencil.C

\*---------------------------------------------------------------------------*/

#ifndef cellStencil_H
#define cellStencil_H

#include "labelList.H"
#include "SubList.H"
#include "vector.H"
#include "Map.H"
#include "globalIndex.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class cellStencil;
Ostream& operator<<(Ostream&, const cellStencil&);
Istream& operator>>(Istream&, cellStencil&);

class cellStencil
:
    public labelList
{
    //- Owner cell (local)
    label owner_;

    //- Current max level
    label maxLevel_;

    //- Location of cell centre
    vector centre_;

    //- Start of each level
    labelList offsets_;

    //- Local cells
    mutable labelList localStencil_;

    //- Local levels
    mutable labelList localOffsets_;


public:

    // Constructors

        //- Construct null
    cellStencil();

        //- Copy constructor
    cellStencil(const cellStencil& stencil);

        //- Construct given inputs
    cellStencil
    (
        const label own,
        const labelList& stencil,
        const labelList& levels,
        const vector& centre
    );

        //- Construct will inputs and level map
    cellStencil
    (
        const label own,
        const labelList& stencil,
        const Map<label>& levelMap,
        const vector& centre
    );

    // Member functions

        // Access

            //- Return if the owner is local
            bool isLocal() const
            {
                return owner_ >= 0;
            }

            //- Return the owner index (global)
            inline label owner() const
            {
                return this->first();
            }

            //- Return the owner index (local)
            inline label localOwner() const
            {
                return owner_;
            }

            //- Return the centre point
            inline const vector& centre() const
            {
               return centre_;
            }

            //- Return the offsets
            inline const labelList& offsets() const
            {
                return offsets_;
            }

            //- Return the full local stencil
            const labelList& localStencil() const
            {
            return
                Pstream::parRun()
              ? localStencil_
              : static_cast<const labelList&>(*this);
            }

            //- Return the local stencil up to the given level
            SubList<label> localStencil(const label level) const;

            //- Return the local stencil between two levels
            SubList<label> localStencil
            (
            const label level0,
                const label level1
            ) const;


    // Edit

        //- Update the local stencil
        void updateLocalStencil(const globalIndex&);


   // Operators

    //- Assignement
    void operator=(const cellStencil& stencil)
    {
        labelList::operator=(stencil);
        owner_ = stencil.owner_;
        maxLevel_ = stencil.maxLevel_;
        centre_ = stencil.centre_;
        offsets_ = stencil.offsets_;
        localStencil_ = stencil.localStencil_;
        localOffsets_ = stencil.localOffsets_;
    }

        //- Return the stencil up to the given level
    SubList<label> operator()(const label level) const;

        //- Return the stencil between two levels
    SubList<label> operator()(const label level0, const label level1) const;

        // Output
    friend Ostream& operator<<(Ostream& os, const cellStencil&);

        // Input
    friend Istream& operator>>(Istream& is, cellStencil&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
