/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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
    Foam::scalarLookupTable3D

Description
    Table used to lookup scalar values given a 3D table

SourceFiles
    scalarLookupTable3D.C

\*---------------------------------------------------------------------------*/

#ifndef tableReader_H
#define tableReader_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fileOperation.H"
#include "stringOps.H"
#include "List2D.H"
#include "List3D.H"
#include "Modifier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef List2D<List<token>> entryTable;
extern HashTable<entryTable> readTables;

//- Remove comments from a line
void removeComments(string& line);

template<class Type>
const dictionary& readComponent
(
    const dictionary& parentDict,
    const word& name,
    autoPtr<Modifier<Type>>& mod,
    List<Type>& values,
    const bool canRead
);

//- Read a 1D table from a file
template<class Type>
void read1DTable
(
    const fileName& file,
    const char delim,
    List<Type>& values,
    const bool determineSize = true
);


//- Read a 2D table
template<class Type>
void read2DTable
(
    const fileName& file,
    const char delim,
    List2D<Type>& data,
    const bool flip,
    const bool determineSize = false
);


const entryTable& read2DTable
(
    const fileName& file,
    const char delim,
    const label startRow = 0,
    const bool flip = false
);

//- Read a column of a table
template<class Type>
List<Type> readColumn
(
    const entryTable& entries,
    const label col
);

//- Read a row of a table
template<class Type>
List<Type> readRow
(
    const entryTable& entries,
    const label col
);

//- Read a 3D table
template<class Type>
void read3DTable
(
    const fileName& file,
    const char delim,
    const char rowDelim,
    List3D<Type>& data,
    const bool flip,
    const bool determineSize = false
);

//- Template specialisation for pTraits<char>
char readChar(Istream&);
Istream& operator>>(Istream&, char&);

template<>
class pTraits<char>
{
    char p_;

public:

    //- Component type
    typedef char cmptType;


    // Static Data Members

        static const char* const typeName;

    // Constructors

        //- Construct from primitive
        explicit pTraits(const char&);

        //- Construct from Istream
        pTraits(Istream&);

    // Member Functions

        //- Access to the char value
        operator char() const
        {
            return p_;
        }

        //- Access to the char value
        operator char&()
        {
            return p_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tableReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
