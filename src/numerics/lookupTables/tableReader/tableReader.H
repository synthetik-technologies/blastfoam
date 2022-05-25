/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technologies
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
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Remove comments from a line
void removeComments(string& line);

template<class Type>
bool readComponent
(
    const dictionary& parentDict,
    const word& name,
    word& modType,
    Field<Type>& values,
    const List<List<string>>& table = List<List<string>>()
);

//- Read a 1D table from a file
template<class Type>
void read1DTable
(
    const fileName& file,
    const string& delim,
    Field<Type>& values,
    const bool determineSize = true
);


//- Read a 2D table
template<class Type, template<class> class ListType>
void read2DTable
(
    const fileName& file,
    const string& delim,
    List<ListType<Type>>& data,
    const bool flip,
    const bool determineSize = false
);


List<List<string>> read2DTable
(
    const fileName& file,
    const string& delim,
    const label startRow = 0,
    const bool flip = false
);

//- Read a column of a table
template<class Type>
List<Type> readColumn
(
    const List<List<string>>& entries,
    const label col
);

//- Read a 3D table
template<class Type, template<class> class ListType1, template<class> class ListType2>
void read3DTable
(
    const fileName& file,
    const string& delim,
    const string& rowDelim,
    List<ListType1<ListType2<Type>>>& data,
    const bool flip,
    const bool determineSize = false
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tableReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
