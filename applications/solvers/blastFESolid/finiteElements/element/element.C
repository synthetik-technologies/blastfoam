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

\*---------------------------------------------------------------------------*/

#include "element.H"
#include "tensor2D.H"
#include "degenerateMatcher.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::element::element()
:
    List<label>(),
    elem_(nullptr),
    ir_(nullptr),
    geoType_(GeoType::UNKNOWN)
{}

Foam::element::element
(
    const primitiveMesh& mesh,
    const label index,
    const label order,
    const label intOrder,
    const GeoType::Type type
)
:
    List<label>(),
    elem_(nullptr),
    ir_(nullptr),
    geoType_(GeoType::UNKNOWN)
{
    initialize(mesh, index, order, intOrder, type);
}


// Foam::element::element(const labelList& elem, const ReadType rt)
//
// :
//     List<label>(),
//     elem_(nullptr),
//     ir_(nullptr),
//     geoType_(GeoType::UNKNOWN)
// {
//     initialize(elem, rt);
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::element::~element()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::element::initialize
(
    const primitiveMesh& mesh,
    const label index,
    const label order,
    const label intOrder,
    const GeoType::Type gt
)
{
    geoType_ = gt;
    switch (geoType_)
    {
        case GeoType::CELL:
        {
            cellShape cs(degenerateMatcher::match(mesh, index));
            elem_ = finiteElement::getRefFiniteElement
            (
                cs.model().name(),
                order
            );
//             switch (elem_->elemType())
//             {
//                 case ElementType::HEX:
//                 {
// //                     faceList faces(cs.model().faces(cs));
// //                     const cell c(identity(faces.size()));
// //
// //                     label bottomFaceI = 0;
// //                     face& bottomFace = faces[bottomFaceI];
// //                     labelList::operator=(bottomFace);
// //
// //                     edgeList edges = c.edges(faces);
// //                     label topStart = -1;
// //                     forAll(edges, ei)
// //                     {
// //                         const edge& e = edges[ei];
// //                         if (e[0] == bottomFace[0])
// //                         {
// //                             if (findIndex(bottomFace, e[1]) < 0)
// //                             {
// //                                 topStart = e[1];
// //                                 break;
// //                             }
// //                         }
// //                         else if (e[1] == bottomFace[0])
// //                         {
// //                             if (findIndex(bottomFace, e[0]) < 0)
// //                             {
// //                                 topStart = e[0];
// //                                 break;
// //                             }
// //                         }
// //                     }
// //
// //                     face topFace(bottomFace.size());
// //                     {
// //                         const label topFaceI =
// //                             c.opposingFaceLabel(bottomFaceI, faces);
// //                         label pi = findIndex(faces[topFaceI], topStart);
// //                         forAll(topFace, i)
// //                         {
// //                             topFace[i] = faces[topFaceI][pi];
// //                             pi = faces[topFaceI].fcIndex(pi);
// //                         }
// //                     }
// //                     labelList::append(topFace);
//                     faceList faces(cs.model().faces(cs));
//                     const cell c(identity(faces.size()));
//
//                     label bottomFaceI = -1;
//                     scalar minz = great;
//                     bool flip = true;
//                     forAll(faces, fi)
//                     {
//                         vector fc(faces[fi].centre(mesh.points()));
//
//                         if (fc.z() < minz)
//                         {
//                             minz = fc.z();
//                             bottomFaceI = fi;
//                         }
//                     }
//                     face& bottomFace = faces[bottomFaceI];
//
//                     {
//                         vector fa(bottomFace.area(mesh.points()));
//                         if (fa.z() < 0) bottomFace.flip();
//
//                         scalarField magPts(mag(pointField(mesh.points(), bottomFace)));
//                         inplaceRotateList
//                         (
//                             bottomFace,
//                             bottomFace.size() - findMin(magPts)
//                         );
//                     }
//                     labelList::operator=(bottomFace);
//
//                     edgeList edges = c.edges(faces);
//                     label topStart = -1;
//                     forAll(edges, ei)
//                     {
//                         const edge& e = edges[ei];
//                         if (e[0] == bottomFace[0])
//                         {
//                             if (findIndex(bottomFace, e[1]) < 0)
//                             {
//                                 topStart = e[1];
//                                 break;
//                             }
//                         }
//                         else if (e[1] == bottomFace[0])
//                         {
//                             if (findIndex(bottomFace, e[0]) < 0)
//                             {
//                                 topStart = e[0];
//                                 break;
//                             }
//                         }
//                     }
//
//                     face topFace(bottomFace.size());
//                     {
//                         const label topFaceI =
//                             c.opposingFaceLabel(bottomFaceI, faces);
//                         label pi = findIndex(faces[topFaceI], topStart);
//                         forAll(topFace, i)
//                         {
//                             topFace[i] = faces[topFaceI][pi];
//                             pi = faces[topFaceI].fcIndex(pi);
//                         }
//
//                         vector fa(topFace.area(mesh.points()));
//                         if (fa.z() < 0) topFace.flip();
//                     }
//                     labelList::append(topFace);
//                     break;
//                 }
//                 default:
                // {
                    this->transfer(cs);
                    // break;
                // }
            // }


            break;
        }
        case GeoType::FACE:
        {
            const face& f = mesh.faces()[index];
            if (f.size() > 4 || f.size() < 3)
            {
                FatalErrorInFunction
                    << "Only triangle or quads are supported for 2D finite elements" << endl
                    << abort(FatalError);
            }
            List<label>::operator=(f);
            elem_ = finiteElement::getRefFiniteElement
            (
                f.size() == 3 ? "tri" : "quad",
                order
            );
            break;
        }
        case GeoType::EDGE:
        {
            const edge& e = mesh.edges()[index];
            List<label>::operator=({e[0], e[1]});
            elem_ = finiteElement::getRefFiniteElement("seg", order);
            break;
        }
        case GeoType::POINT:
        {
            this->setSize(1);
            List<label>::operator=(index);
            elem_ = finiteElement::getRefFiniteElement("pt", order);
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Geometric type has not been set" << endl
                << abort(FatalError);
        }
    }
//     ir_ = &elem_->ir();
    ir_ = &integrationRule::getRule(type(), intOrder);
}


// void Foam::element::initialize(const labelList& elem, const ReadType rt);
// {
// }


void Foam::element::initializeShell
(
    const primitiveMesh& mesh,
    const label index,
    const label order,
    const label orderRS,
    const label orderT,
    const GeoType::Type gt
)
{
    geoType_ = gt;
    switch (geoType_)
    {
        case GeoType::FACE:
        {
            const face& f = mesh.faces()[index];
            if (f.size() > 4 || f.size() < 3)
            {
                FatalErrorInFunction
                    << "Only triangle or quads are supported for 2D finite elements" << endl
                    << abort(FatalError);
            }
            List<label>::operator=(f);
            elem_ = shellFiniteElement::getRefFiniteElement
            (
                f.size() == 3 ? "tri" : "quad",
                order
            );
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Geometric type has not been set" << endl
                << abort(FatalError);
        }
    }
//     ir_ = &elem_->ir();
    ir_ = &shellIntegrationRule::getRule(type(), orderRS, orderT);
}

Foam::scalarList Foam::element::calcShape(const vector& pt) const
{
    return elem_->calcShape(pt);
}


Foam::scalarRectangularMatrix Foam::element::calcDShape
(
    const vector& pt
) const
{
    return elem_->calcDShape(pt);
}


Foam::tensor Foam::element::calcJ
(
    const vector& ip,
    const pointField& meshNodes,
    const scalarRectangularMatrix& dshape
) const
{
    UIndirectList<vector> nodes(meshNodes, *this);
    tensor J(Zero);
    forAll(nodes, nodei)
    {
        for (label dimi = 0; dimi < 3; dimi++)
        {
            for (label dimj = 0; dimj < dshape.n(); dimj++)
            {
                J(dimi, dimj) += nodes[nodei][dimi]*dshape(nodei, dimj);
            }
        }
    }
    return J;
}


Foam::tensor Foam::element::calcJ
(
    const vector& ip,
    const pointField& meshNodes
) const
{
    return calcJ(ip, meshNodes, calcDShape(ip));
}


Foam::tensor Foam::element::calcInvJ
(
    const vector& ip,
    const pointField& meshNodes
) const
{
    if (geoType_ == GeoType::POINT)
    {
        return tensor::zero;
    }

    tensor J(calcJ(ip, meshNodes));
    return calcInvJ(J);
}


Foam::tensor Foam::element::calcInvJ
(
    const tensor& J
) const
{
    if (geoType_ == GeoType::POINT)
    {
        return tensor::zero;
    }
    else if (geoType_ == GeoType::CELL)
    {
        return J.inv();
    }

    label valid[3] = {-1, -1, -1};
    label vi = 0;
    for (label dimj = 0; dimj < 3; dimj++)
    {
        bool allZero = true;
        for (label dimi = 0; dimi < 3; dimi++)
        {
            if (mag(J(dimi, dimj)) > small)
            {
                allZero = false;
                break;
            }
        }
        if (!allZero)
        {
            valid[vi++] = dimj;
        }
    }

    if (vi == 3)
    {
        return J.inv();
    }
    else if (vi == 2)
    {
        const label i = valid[0];
        const label j = valid[1];
        scalar E = sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i));
        scalar G = sqr(J(0, j)) + sqr(J(1, j)) + sqr(J(2, j));
        scalar F =
            J(0, i)*J(0, j)
          + J(1, i)*J(1, j)
          + J(2, i)*J(2, j);
        const scalar T = 1.0/E*G - F*F;
        E *= T;
        G *= T;
        F *= T;

        tensor invJ(Zero);
        invJ(i, 0) = J(0, i)*G - J(0, j)*F;
        invJ(j, 0) = J(0, j)*E - J(0, i)*F;
        invJ(i, 1) = J(1, i)*G - J(1, j)*F;
        invJ(j, 1) = J(1, j)*E - J(1, i)*F;
        invJ(i, 2) = J(2, i)*G - J(2, j)*F;
        invJ(j, 2) = J(2, j)*E - J(2, i)*F;
        return invJ;
    }
    else if (vi == 1)
    {
        const label i = valid[0];
        const scalar T = 1.0/sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i));
        tensor invJ(Zero);
        invJ(i, 0) = J(0, i)*T;
        invJ(i, 1) = J(1, i)*T;
        invJ(i, 2) = J(2, i)*T;
        return invJ;
    }
    return tensor::zero;
}


Foam::scalar Foam::element::calcW
(
    const vector& ip,
    const pointField& meshNodes
) const
{
    if (geoType_ == GeoType::POINT)
    {
        return 0.0;
    }

    tensor J(calcJ(ip, meshNodes));
    return calcW(J);
}


Foam::scalar Foam::element::calcW
(
    const tensor& J
) const
{
    if (geoType_ == GeoType::POINT)
    {
        return 0.0;
    }
    else if (geoType_ == GeoType::CELL)
    {
        return det(J);
    }

    label valid[3] = {-1, -1, -1};
    label vi = 0;
    for (label dimj = 0; dimj < 3; dimj++)
    {
        bool allZero = true;
        for (label dimi = 0; dimi < 3; dimi++)
        {
            if (mag(J(dimi, dimj)) > small)
            {
                allZero = false;
                break;
            }
        }
        if (!allZero)
        {
            valid[vi++] = dimj;
        }
    }

    if (vi == 3)
    {
        return det(J);
    }
    else if (vi == 2)
    {
        const label i = valid[0];
        const label j = valid[1];
        const scalar E = sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i));
        const scalar G = sqr(J(0, j)) + sqr(J(1, j)) + sqr(J(2, j));
        const scalar F =
            J(0, i)*J(0, j)
          + J(1, i)*J(1, j)
          + J(2, i)*J(2, j);
        return sqrt(E*G - F*F);
    }
    else if (vi == 1)
    {
        const label i = valid[0];
        return sqrt(sqr(J(0, i)) + sqr(J(1, i)) + sqr(J(2, i)));
    }
    return 0.0;
}


Foam::vector Foam::element::calcOrtho(const tensor& J) const
{
//     label valid[3] = {-1, -1, -1};
//     label vi = 0;
//     for (label dimj = 0; dimj < 3; dimj++)
//     {
//         bool allZero = true;
//         for (label dimi = 0; dimi < 3; dimi++)
//         {
//             if (mag(J(dimi, dimj)) > small)
//             {
//                 allZero = false;
//                 break;
//             }
//         }
//         if (!allZero)
//         {
//             valid[vi++] = dimj;
//         }
//     }
//
//     if (vi == 3)
    {
        return
            normalised
            (
//                 vector
//                 (
//                     J[1]*J[5] - J[2]*J[4],
//                     J[2]*J[3] - J[0]*J[5],
//                     J[0]*J[4] - J[1]*J[3]
//                 )
                vector
                (
                    J(1, 0)*J(2, 1) - J(2, 0)*J(1, 1),
                    J(2, 0)*J(0, 1) - J(0, 0)*J(2, 1),
                    J(0, 0)*J(1, 1) - J(1, 0)*J(0, 1)
                )
            );
    }
//     else if (vi == 2)
//     {
//         const label i = valid[0];
//         const label j = valid[1];
//         vector n(Zero);
//         n[i] = -J(i, j);
//         n[j] = J(i, i);
//         return normalised(n);
//     }
//     else if (vi == 1)
//     {
//         NotImplemented;
//     }
//     return vector::zero;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //


// ************************************************************************* //
