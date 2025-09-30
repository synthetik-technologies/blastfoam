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

#include "segmentFiniteElement.H"
#include "GaussianQuadrature.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(segment, 0);
    addToRunTimeSelectionTable(finiteElement, segment, type);
    addNamedToRunTimeSelectionTable(finiteElement, segment, type, seg);

    addFE(seg2); addFEMap(seg2, msh, 1);
    addFE(seg3); addFEMap(seg3, msh, 8);
    addFE(seg4); addFEMap(seg4, msh, 26);
    addFE(seg5); addFEMap(seg5, msh, 27);
    addFE(seg6); addFEMap(seg6, msh, 28);
    addFE(seg7); addFEMap(seg7, msh, 62);
    addFE(seg8); addFEMap(seg8, msh, 63);
    addFE(seg9); addFEMap(seg9, msh, 64);
    addFE(seg10); addFEMap(seg10, msh, 65);
    addFE(seg11); addFEMap(seg11, msh, 66);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::segment::segment(const label order)
:
    FiniteElementBase<ElementType::SEG>(order, order+1),
    w_(this->nNodes(), 0.0)
{
    CmptList<vector> x(nodes_, vector::X);
    if (order == 0)
    {
        x[0] = 0.5;
        w_ = 1.0;
    }
    else
    {
//         List<scalar> _x;
//         GaussianQuadrature::calcLobatto(order+1, _x, w_);
//         forAll(x, i)
//         {
//             x[i] = _x[i];
//         }

        scalar dx = 1.0/scalar(order);
        forAll(x, i)
        {
            x[i] = dx*scalar(i);
        }
        w_ = 1.0;
        for (label i = 0; i <= order; i++)
        {
            for (label j = 0; j < i; j++)
            {
                scalar xij = x[i] - x[j];
                w_[i] *= xij;
                w_[j] *= -xij;
            }
        }
        forAll(w_, i)
        {
            w_[i] = 1.0/w_[i];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::segment::~segment()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::segment::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize(cp.size(), vector::zero);
    G_.setSize(order_);
    E_.setSize(0);

    forAll(cp, i)
    {
        P_[i].x() = cp[i];
    }
    for (label i = 0; i <= order_; i++)
    {
        G_[i] = {i, i+1};
    }
}


Foam::tmp<Foam::pointField> Foam::finiteElements::segment::getNodes
(
    const List<vector>& verts
) const
{
    tmp<pointField> tnodes(new pointField(this->nNodes()));
    pointField& nodes = tnodes.ref();

    const CmptList<vector> x(nodes_, vector::X);
    const vector& p0(verts[0]);
    const vector dx(verts[1] - verts[0]);
    for (label i = 0; i < x.size(); i++)
    {
        nodes[i] = p0 + x[i]*dx;
    }
    return tnodes;
}


void Foam::finiteElements::segment::calcShape
(
    const scalar x,
    scalar* shape
) const
{
    if (order_ == 0)
    {
        shape[0] = 1.0;
        return;
    }

    const CmptList<vector> xs(nodes_, vector::X);
    const label p = xs.size()-1;
    scalar sk = 1.0;
    label k, i;
    for (k = 0; k < p; k++)
    {
        if (x >= (xs[k] + xs[k+1])*0.5)
        {
            sk *= x - xs[k];
        }
        else
        {
            for (i = k+1; i <= p; i++)
            {
                sk *= x - xs[i];
            }
            break;
        }
    }
    scalar s = sk*(x - xs[k]);

    for (i = 0; i < k; i++)
    {
        shape[i] = s*w_[i]/(x - xs[i]);
    }
    shape[k] = sk*w_[k];
    for (i++; i <= p; i++)
    {
        shape[i] = s*w_[i]/(x - xs[i]);
    }
}


Foam::scalarList Foam::finiteElements::segment::calcShape
(
    const scalar x
) const
{
    scalarList shape(order_ + 1);
    calcShape(x, shape.data());
    return shape;
}

Foam::scalarList Foam::finiteElements::segment::calcShape
(
    const vector& pt
) const
{
    scalarList shape(order_ + 1);
    calcShape(pt.x(), shape.data());
    return shape;
}


void Foam::finiteElements::segment::calcDShape
(
    const scalar x,
    scalar* dshape
) const
{
    if (order_ == 0)
    {
        dshape[0] = 0.0;
        return;
    }

    List<scalar> shape(order_+1, 0.0);

    const CmptList<vector> xs(nodes_, vector::X);
    scalar sk = 1.0;
    label k, i;
    for (k = 0; k < order_; k++)
    {
        if (x >= (xs[k] + xs[k+1])*0.5)
        {
            sk *= x - xs[k];
        }
        else
        {
            for (i = k+1; i <= order_; i++)
            {
                sk *= x - xs[i];
            }
            break;
        }
    }
    scalar s = sk*(x - xs[k]);

    scalar di = 0.0;
    scalar dk = 0.0;
    for (i = 0; i < k; i++)
    {
        di = 1.0/(x - xs[i]);
        dk += di;
        shape[i] = s*w_[i]*di;
    }
    shape[k] = sk*w_[k];
    for (i++; i <= order_; i++)
    {
        di = 1.0/(x - xs[i]);
        dk += di;
        shape[i] = s*w_[i]*di;
    }
    scalar sp = s*dk + sk;

    for (i = 0; i < k; i++)
    {
        dshape[i] = (sp*w_[i] - shape[i])/(x - xs[i]);
    }
    dshape[k] = dk*shape[k];
    for (i++; i <= order_; i++)
    {
        dshape[i] = (sp*w_[i] - shape[i])/(x - xs[i]);
    }
}


Foam::scalarList Foam::finiteElements::segment::calcDShape
(
    const scalar x
) const
{
    scalarList dshape(order_+1);
    calcDShape(x, dshape.data());
    return dshape;
}


Foam::scalarRectangularMatrix Foam::finiteElements::segment::calcDShape
(
    const vector& pt
) const
{
    scalarRectangularMatrix dshape(order_+1, this->nDims());
    calcDShape(pt.x(), dshape.v());
    return dshape;
}


void Foam::finiteElements::segment::calcDShape
(
    const scalar x,
    scalarList& shape,
    scalarList& dshape
) const
{
    shape.setSize(this->nNodes());
    dshape.setSize(this->nNodes(), 1);
    if (order_ == 0)
    {
        shape = 1.0;
        dshape = 0.0;
        return;
    }

    const CmptList<vector> xs(nodes_, vector::X);
    scalar sk = 1.0;
    label k, i;
    for (k = 0; k < order_; k++)
    {
        if (x >= (xs[k] + xs[k+1])*0.5)
        {
            sk *= x - xs[k];
        }
        else
        {
            for (i = k+1; i <= order_; i++)
            {
                sk *= x - xs[i];
            }
            break;
        }
    }
    scalar s = sk*(x - xs[k]);

    scalar di = 0.0;
    scalar dk = 0.0;
    for (i = 0; i < k; i++)
    {
        di = 1.0/(x - xs[i]);
        dk += di;
        shape[i] = s*w_[i]*di;
    }
    shape[k] = sk*w_[k];
    for (i++; i <= order_; i++)
    {
        di = 1.0/(x - xs[i]);
        dk += di;
        shape[i] = s*w_[i]*di;
    }
    scalar sp = s*dk + sk;

    for (i = 0; i < k; i++)
    {
        dshape[i] = (sp*w_[i] - shape[i])/(x - xs[i]);
    }
    dshape[k] = dk*shape[k];
    for (i++; i <= order_; i++)
    {
        dshape[i] = (sp*w_[i] - shape[i])/(x - xs[i]);
    }
}


Foam::label Foam::finiteElements::segment::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 3; // SEGMENT
        case 2:
            return 21; // QUADRATIC_SEGMENT
        default:
            return 68; // LAGRANGE_SEGMENT
    }
}

Foam::label Foam::finiteElements::segment::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 1;
        case 2:
            return 8;
        default:
            return -1;
    }
}


void Foam::finiteElements::segment::vtkData
(
    labelList& data,
    label& start,
    const labelList& labels
) const
{
    switch (order_)
    {
//         case 0:
//         case 1:
//         {
//             forAll(labels, i)
//             {
//                 data[start++] = labels[i];
//             }
//             return;
//         }
        default:
        {
            forAll(labels, i)
            {
                data[start++] = labels[i];
            }
            return;
        }
    }
}


// ************************************************************************* //
