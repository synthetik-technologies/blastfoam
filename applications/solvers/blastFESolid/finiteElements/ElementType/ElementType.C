/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024
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

\*---------------------------------------------------------------------------*/

#include "ElementType.H"
#include "triFace.H"
#include "cell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
template<>
const char* Foam::NamedEnum<Foam::ElementType::Type, 8>::names[] =
{
    "pt",

    "seg",

    "tri",
    "quad",

    "tet",
    "hex",
    "prism",
    "pyr"
};

const Foam::NamedEnum<Foam::ElementType::Type, 8> Foam::ElementType::cellModelNames;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

//- Point specializations
Foam::edgeList Foam::ElementType::edges
(
    const Type t,
    const List<label>& elem
)
{
    switch (t)
    {
        case PT:
            return Geometry<PT>::edges(elem);
        case SEG:
            return Geometry<SEG>::edges(elem);
        case TRI:
            return Geometry<TRI>::edges(elem);
        case QUAD:
            return Geometry<QUAD>::edges(elem);
        case TET:
            return Geometry<TET>::edges(elem);
        case HEX:
            return Geometry<HEX>::edges(elem);
        case PRISM:
            return Geometry<PRISM>::edges(elem);
        case PYR:
            return Geometry<PYR>::edges(elem);
        default:
            return edgeList();
    }
}


Foam::faceList Foam::ElementType::faces
(
    const Type t,
    const List<label>& elem
)
{
    switch (t)
    {
        case PT:
            return Geometry<PT>::faces(elem);
        case SEG:
            return Geometry<SEG>::faces(elem);
        case TRI:
            return Geometry<TRI>::faces(elem);
        case QUAD:
            return Geometry<QUAD>::faces(elem);
        case TET:
            return Geometry<TET>::faces(elem);
        case HEX:
            return Geometry<HEX>::faces(elem);
        case PRISM:
            return Geometry<PRISM>::faces(elem);
        case PYR:
            return Geometry<PYR>::faces(elem);
        default:
            return faceList();
    }
}


// void Foam::ElementType::toRef
// (
//     const Type t,
//     labelList& elem,
//     const faceList& fs,
//     const pointField& pts
// )
// {}


//- Point specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::PT>::edges
(
    const List<label>& elem
)
{
    return edgeList();
}


Foam::faceList Foam::Geometry<Foam::ElementType::PT>::faces
(
    const List<label>& elem
)
{
    return faceList();
}


void Foam::Geometry<Foam::ElementType::PT>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}


//- Segment specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::SEG>::edges
(
    const List<label>& elem
)
{
    return edgeList({edge(elem[0], elem[1])});
}


Foam::faceList Foam::Geometry<Foam::ElementType::SEG>::faces
(
    const List<label>& elem
)
{
    return faceList();
}


void Foam::Geometry<Foam::ElementType::SEG>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}


//- Triangle specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::TRI>::edges
(
    const List<label>& elem
)
{
    return edgeList
    (
        {
            edge(elem[0], elem[1]),
            edge(elem[1], elem[2]),
            edge(elem[2], elem[0])
        }
    );
}


Foam::faceList Foam::Geometry<Foam::ElementType::TRI>::faces
(
    const List<label>& elem
)
{
    return faceList({triFace({elem[0], elem[1], elem[2]})});
}


void Foam::Geometry<Foam::ElementType::TRI>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}


//- Quadrilateral specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::QUAD>::edges
(
    const List<label>& elem
)
{
    return edgeList
    (
        {
            edge(elem[0], elem[1]),
            edge(elem[1], elem[2]),
            edge(elem[2], elem[3]),
            edge(elem[3], elem[0])
        }
    );
}


Foam::faceList Foam::Geometry<Foam::ElementType::QUAD>::faces
(
    const List<label>& elem
)
{
    return faceList({face({elem[0], elem[1], elem[2], elem[3]})});
}


void Foam::Geometry<Foam::ElementType::QUAD>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}


//- Tetrahedron specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::TET>::edges
(
    const List<label>& elem
)
{
    return edgeList
    (
        {
            edge(elem[0], elem[1]),
            edge(elem[0], elem[2]),
            edge(elem[0], elem[3]),
            edge(elem[1], elem[2]),
            edge(elem[1], elem[3]),
            edge(elem[2], elem[3])
        }
    );
}


Foam::faceList Foam::Geometry<Foam::ElementType::TET>::faces
(
    const List<label>& elem
)
{
    return faceList
    (
        {
            triFace({elem[1], elem[2], elem[3]}),
            triFace({elem[0], elem[3], elem[2]}),
            triFace({elem[0], elem[1], elem[3]}),
            triFace({elem[0], elem[2], elem[1]})
        }
    );
}


void Foam::Geometry<Foam::ElementType::TET>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}


//- Hexahedron specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::HEX>::edges
(
    const List<label>& elem
)
{
    return edgeList
    (
        {
            edge(elem[0], elem[1]),
            edge(elem[1], elem[2]),
            edge(elem[3], elem[2]),
            edge(elem[0], elem[3]),
            edge(elem[4], elem[5]),
            edge(elem[5], elem[6]),
            edge(elem[7], elem[6]),
            edge(elem[4], elem[7]),
            edge(elem[0], elem[4]),
            edge(elem[1], elem[5]),
            edge(elem[2], elem[6]),
            edge(elem[3], elem[7])
        }
    );
}


Foam::faceList Foam::Geometry<Foam::ElementType::HEX>::faces
(
    const List<label>& elem
)
{
    return faceList
    (
        {
            face({elem[3], elem[2], elem[1], elem[0]}),
            face({elem[0], elem[1], elem[5], elem[4]}),
            face({elem[1], elem[2], elem[6], elem[5]}),
            face({elem[2], elem[3], elem[7], elem[6]}),
            face({elem[3], elem[0], elem[4], elem[7]}),
            face({elem[4], elem[5], elem[6], elem[7]})
        }
    );
}


void Foam::Geometry<Foam::ElementType::HEX>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{
    const cell c(identity(fs.size()));

    label bottomFaceI = -1;
    scalar minz = great;
    forAll(fs, fi)
    {
        vector fc(fs[fi].centre(pts));

        if (fc.z() < minz)
        {
            minz = fc.z();
            bottomFaceI = fi;
        }
    }
    face bottomFace(fs[bottomFaceI]);

    {
        vector fa(bottomFace.area(pts));
        if (fa.z() < 0) bottomFace.flip();

        scalarField magPts(mag(pointField(pts, bottomFace)));
        inplaceRotateList<List, label>
        (
            bottomFace,
            bottomFace.size() - findMin(magPts)
        );
    }

    elem = bottomFace;

    edgeList edges = c.edges(fs);
    label topStart = -1;
    forAll(edges, ei)
    {
        const edge& e = edges[ei];
        if (e[0] == bottomFace[0])
        {
            if (findIndex(bottomFace, e[1]) < 0)
            {
                topStart = e[1];
                break;
            }
        }
        else if (e[1] == bottomFace[0])
        {
            if (findIndex(bottomFace, e[0]) < 0)
            {
                topStart = e[0];
                break;
            }
        }
    }

    face topFace(bottomFace.size());
    {
        const label topFaceI =
            c.opposingFaceLabel(bottomFaceI, fs);
        label pi = findIndex(fs[topFaceI], topStart);
        forAll(topFace, i)
        {
            topFace[i] = fs[topFaceI][pi];
            pi = fs[topFaceI].fcIndex(pi);
        }

        vector fa(topFace.area(pts));
        if (fa.z() < 0) topFace.flip();
    }
    elem.append(topFace);
}

//- Prism specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::PRISM>::edges
(
    const List<label>& elem
)
{
    return edgeList
    (
        {
            edge(elem[0], elem[1]),
            edge(elem[1], elem[2]),
            edge(elem[2], elem[0]),
            edge(elem[3], elem[4]),
            edge(elem[4], elem[5]),
            edge(elem[5], elem[3]),
            edge(elem[0], elem[3]),
            edge(elem[1], elem[4]),
            edge(elem[2], elem[5])
        }
    );
}


Foam::faceList Foam::Geometry<Foam::ElementType::PRISM>::faces
(
    const List<label>& elem
)
{
    return faceList
    (
        {
            triFace({elem[3], elem[2], elem[1]}),
            triFace({elem[3], elem[4], elem[5]}),
            face({elem[0], elem[1], elem[4], elem[3]}),
            face({elem[1], elem[2], elem[5], elem[4]}),
            face({elem[2], elem[0], elem[3], elem[5]})
        }
    );
}


void Foam::Geometry<Foam::ElementType::PRISM>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}


//- Pyramid specializations
Foam::edgeList Foam::Geometry<Foam::ElementType::PYR>::edges
(
    const List<label>& elem
)
{
    return edgeList
    (
        {
            edge(elem[0], elem[1]),
            edge(elem[1], elem[2]),
            edge(elem[2], elem[3]),
            edge(elem[3], elem[0]),
            edge(elem[0], elem[4]),
            edge(elem[1], elem[4]),
            edge(elem[2], elem[4]),
            edge(elem[3], elem[4])
        }
    );
}


Foam::faceList Foam::Geometry<Foam::ElementType::PYR>::faces
(
    const List<label>& elem
)
{
    return faceList
    (
        {
            face({elem[3], elem[2], elem[1], elem[0]}),
            triFace({elem[0], elem[1], elem[4]}),
            triFace({elem[1], elem[2], elem[4]}),
            triFace({elem[2], elem[3], elem[4]}),
            triFace({elem[3], elem[0], elem[4]})
        }
    );
}


void Foam::Geometry<Foam::ElementType::PYR>::toRef
(
    labelList& elem,
    const faceList& fs,
    const pointField& pts
)
{}

// Foam::label Foam::ElementType::nVertices(const Type type)
// {
//     switch (type)
//     {
//         case PT: return 1;
//         case SEG: return 2;
//         case TRI: return 3;
//         case QUAD: return 4;
//         case TET: return 4;
//         case HEX: return 8;
//         case PRISM: return 6;
//         case PYR: return 5;
//         default:
//             NotImplemented;
//             return -1;
//     }
// }
//
// Foam::label Foam::ElementType::nEdges(const Type type)
// {
//     switch (type)
//     {
//         case PT: return 0;
//         case SEG: return 1;
//         case TRI: return 3;
//         case QUAD: return 4;
//         case TET: return 6;
//         case HEX: return 12;
//         case PRISM: return 9;
//         case PYR: return 8;
//         default: return -1;
//     }
// }
//
// Foam::label Foam::ElementType::nFaces(const Type type)
// {
//     switch (type)
//     {
//         case PT:
//         case SEG: return 0;
//         case TRI:
//         case QUAD: return 1;
//         case TET: return 4;
//         case HEX: return 6;
//         case PRISM: return 5;
//         case PYR: return 5;
//         default: return -1;
//     }
// }
//
// Foam::label Foam::ElementType::nNodes
// (
//     const Type type,
//     const label order
// )
// {
//     switch (type)
//     {
//         case PT: return 1;
//         case SEG: return order+1;
//         case TRI: return (order+1)*(order+2)/2;
//         case QUAD: return (order+1)*(order+1);
//         case TET:
//             return (order+1)*(order+2)*(order+3)/6;
//         case HEX:
//             return (order+1)*(order+1)*(order+1);
//         case PRISM:
//             return (order+1)*(order+1)*(order+2)/2;
//         case PYR:
//             return 4; // only linear
//         default: return -1;
//     }
// }
//
// Foam::label Foam::ElementType::nEdgeNodes(const Type type, const label order)
// {
//     switch (type)
//     {
//         case PT: return 0;
//         default: return order - 1;
//     }
// }
//
// Foam::label Foam::ElementType::nFaceNodes
// (
//     const Type type,
//     const label order,
//     const label fi
// )
// {
//     switch (type)
//     {
//         case PT:
//         case SEG: return 0;
//         case QUAD: return (order-1)*(order-1);
//         case TRI:
//         case TET:
//         {
//             label n = 0;
//             for (label i = 2; i < order; i++) n += order - i;
//             return n;
//         }
//         case HEX: return (order-1)*(order-1);
//         case PRISM:
//             return
//                 fi < 2
//                 ? nFaceNodes(TRI, order, 0)
//                 : (order-1)*(order-1);
//         case PYR:
//             return
//                 fi == 0
//                 ? (order-1)*(order-1)
//                 : nFaceNodes(TRI, order, 0);
//         default: return -1;
//     }
// }
//
// Foam::label Foam::ElementType::nInternalNodes
// (
//     const Type type,
//     const label order
// )
// {
//     switch (type)
//     {
//         case PT:
//         case SEG:
//         case TRI:
//         case QUAD: return 0;
//         case TET:
//         {
//             label n = 0;
//             for (label k = 1; k < order; k++)
//                 for (label j = 1; j+k < order; j++)
//                     for (label i = 1; i+j+k < order; i++)
//                         n++;
//             return n;
//         }
//         case HEX:
//             return (order-1)*(order-1)*(order-1);
//         case PRISM:
//             return nFaceNodes(TRI, order, 0)*(order-1);
//         case PYR:
//             return 0; // only linear
//         default: return -1;
//     }
// }
