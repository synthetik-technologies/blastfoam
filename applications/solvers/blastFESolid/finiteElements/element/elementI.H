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

#include "element.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::GeoType::Type Foam::element::geoType() const
{
    return geoType_;
}


inline Foam::ElementType::Type Foam::element::type() const
{
    return elem_->elemType();
}


inline Foam::label Foam::element::nDims() const
{
    return elem_->nDims();
}


inline Foam::label Foam::element::order() const
{
    return elem_->order();
}


inline const Foam::integrationRule& Foam::element::ir() const
{
    return *ir_;
}


inline Foam::label Foam::element::nVertices() const
{
    return elem_->nVertices();
}


inline Foam::label Foam::element::nEdges() const
{
    return elem_->nEdges();
}


inline Foam::label Foam::element::nFaces() const
{
    return elem_->nFaces();
}


inline Foam::label Foam::element::nNodes() const
{
    return elem_->nNodes();
}


inline Foam::label Foam::element::nEdgeNodes() const

{
    return elem_->nEdgeNodes();
}


inline Foam::label Foam::element::nFaceNodes(const label fi) const
{
    return elem_->nFaceNodes(fi);
}


inline Foam::label Foam::element::nInternalNodes() const
{
    return elem_->nInternalNodes();
}


inline const Foam::pointField& Foam::element::nodes() const
{
    return elem_->nodes();
}


inline const Foam::finiteElement& Foam::element::elem() const
{
    return *elem_;
}

inline Foam::edgeList Foam::element::edges() const
{
    return ElementType::edges(type(), *this);
}


inline Foam::faceList Foam::element::faces() const
{
    return ElementType::faces(type(), *this);
}


inline Foam::tmp<Foam::pointField> Foam::element::getNodes
(
    const List<vector>& verts
) const
{
    return elem_->getNodes(verts);
}


inline void Foam::element::setIntegrationRule(const integrationRule& ir)
{
    ir_ = &ir;
}


inline void Foam::element::setIntegrationOrder(const label order)
{
    ir_ = &integrationRule::getRule
    (
        type(),
        order
    );
}

inline void Foam::element::setFiniteElement(const finiteElement& fe)
{
    elem_ = &fe;
}


inline void Foam::element::setGeoType(const GeoType::Type& gt)
{
    geoType_ = gt;
}

// ************************************************************************* //
