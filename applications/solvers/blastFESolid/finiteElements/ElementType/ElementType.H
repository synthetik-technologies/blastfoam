/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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
    Foam::ElementType

Description
    Types of elements and corresponding geometric information

SourceFiles
    ElementType.C

\*---------------------------------------------------------------------------*/
#ifndef ElementType_H
#define ElementType_H

#include "NamedEnum.H"
#include "edge.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class ElementType Declaration
\*---------------------------------------------------------------------------*/

class ElementType
{
public:

    //- Enumeration of element types
    enum Type
    {
        UNKNOWN = -1,

        PT = 0,

        SEG = 1,

        TRI = 2,
        QUAD = 3,

        TET  = 4,
        HEX = 5,
        PRISM = 6,
        PYR = 7,

        SIZE = 8
    };

    //- Names corresponding to cellModels
    static const NamedEnum<Type, 8> cellModelNames;

    //- Return edges
    static edgeList edges(const Type t, const List<label>& elem);

    //- Return faces
    static faceList faces(const Type t, const List<label>& elem);

//     //- Convert shape to reference (ordering by coordinates)
//     void toRef
//     (
//         const Type,
//         labelList&,
//         const faceList&,
//         const pointField&
//     );
};


template<ElementType::Type T>
class Geometry
{
public:

    //- Number of dimensions
    static const label nDims;

    //- Number of vertices
    static const label nVertices;

    //- Number of edges
    static const label nEdges;

    //- Number of faces
    static const label nFaces;
};


template<>
class Geometry<ElementType::PT>
{
public:

    //- Number of dimensions
    static const label nDims = 0;

    //- Number of vertices
    static const label nVertices = 1;

    //- Number of edges
    static const label nEdges = 0;

    //- Number of faces
    static const label nFaces = 0;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::SEG>
{
public:

    //- Number of dimensions
    static const label nDims = 1;

    //- Number of vertices
    static const label nVertices = 2;

    //- Number of edges
    static const label nEdges = 1;

    //- Number of faces
    static const label nFaces = 0;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::TRI>
{
public:

    //- Number of dimensions
    static const label nDims = 2;

    //- Number of vertices
    static const label nVertices = 3;

    //- Number of edges
    static const label nEdges = 3;

    //- Number of faces
    static const label nFaces = 1;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::QUAD>
{
public:

    //- Number of dimensions
    static const label nDims = 2;

    //- Number of vertices
    static const label nVertices = 4;

    //- Number of edges
    static const label nEdges = 4;

    //- Number of faces
    static const label nFaces = 1;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::TET>
{
public:

    //- Number of dimensions
    static const label nDims = 3;

    //- Number of vertices
    static const label nVertices = 4;

    //- Number of edges
    static const label nEdges = 6;

    //- Number of faces
    static const label nFaces = 4;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::HEX>
{
public:

    //- Number of dimensions
    static const label nDims = 3;

    //- Number of vertices
    static const label nVertices = 8;

    //- Number of edges
    static const label nEdges = 12;

    //- Number of faces
    static const label nFaces = 6;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::PRISM>
{
public:

    //- Number of dimensions
    static const label nDims = 3;

    //- Number of vertices
    static const label nVertices = 6;

    //- Number of edges
    static const label nEdges = 9;

    //- Number of faces
    static const label nFaces = 5;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


template<>
class Geometry<ElementType::PYR>
{
public:

    //- Number of dimensions
    static const label nDims = 3;

    //- Number of vertices
    static const label nVertices = 5;

    //- Number of edges
    static const label nEdges = 8;

    //- Number of faces
    static const label nFaces = 5;

    //- Return edges
    static edgeList edges(const List<label>& elem);

    //- Return faces
    static faceList faces(const List<label>& elem);

    //- Convert shape to reference (ordering by coordinates)
    void toRef(labelList&, const faceList&, const pointField&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
