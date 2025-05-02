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

#include "triangleFiniteElement.H"
#include "CmptList.H"
#include "LUscalarMatrix.H"
#include "GaussianQuadrature.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(triangle, 0);

    addToRunTimeSelectionTable(finiteElement, triangle, type);
    addNamedToRunTimeSelectionTable(finiteElement, triangle, type, tri);

    typedef ShellFiniteElement<triangle> triangleShell;
    defineTemplateTypeNameAndDebug(triangleShell, 0);

    addToRunTimeSelectionTable(shellFiniteElement, triangleShell, type);
    addNamedToRunTimeSelectionTable(shellFiniteElement, triangleShell, type, tri);

    addFE(tri3); addFEMap(tri3, msh, 2); addShellFE(tri3);
    addFE(tri6); addFEMap(tri6, msh, 9); addShellFE(tri6);
    addFE(tri10); addFEMap(tri10, msh, 21); addShellFE(tri10);
    addFE(tri15); addFEMap(tri15, msh, 23); addShellFE(tri15);
    addFE(tri21); addFEMap(tri21, msh, 25); addShellFE(tri21);
    addFE(tri28); addFEMap(tri28, msh, 42); addShellFE(tri28);
    addFE(tri36); addFEMap(tri36, msh, 43); addShellFE(tri36);
    addFE(tri45); addFEMap(tri45, msh, 44); addShellFE(tri45);
    addFE(tri55); addFEMap(tri55, msh, 45); addShellFE(tri55);
    addFE(tri66); addFEMap(tri66, msh, 46); addShellFE(tri66);
}
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::label Foam::finiteElements::triangle::calcNFaceNodes(const label o)
{
    label n = 0;
    for (label i = 2; i < o; i++) n += o - i;
    return n;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::triangle::triangle(const label order)
:
    FiniteElement<ElementType::TRI>
    (
        order,
        (order + 1)*(order + 2)/2
    ),
    seg_
    (
        dynamicCast<const segment>
        (
            *finiteElement::getRefFiniteElement(ElementType::SEG, order))
    ),
    nFaceNodes_(calcNFaceNodes(order_)),
    invT_()
{
    label ni = 0;

    const CmptList<vector> x(seg_.nodes(), vector::X);
    nodes_[0][0] = x[0];
    nodes_[0][1] = x[0];

    nodes_[1][0] = x[order];
    nodes_[1][1] = x[0];

    nodes_[2][0] = x[0];
    nodes_[2][1] = x[order];

    ni = 3;
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[i];
        nodes_[ni][1] = x[0];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[order-1];
        nodes_[ni][1] = x[i];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[0];
        nodes_[ni][1] = x[order-1];
        ni++;
    }

    for (label j = 1; j < order; j++)
    {
        for (label i = 1; i+j < order; i++)
        {
            scalar w = x[i] + x[j] + x[order-i-j];
            nodes_[ni][0] = x[i]/w;
            nodes_[ni][1] = x[j]/w;
            ni++;
        }
    }

    scalarSquareMatrix T(this->nNodes(), 0.0);
    invT_.setSize(this->nNodes());

    List<scalar> shape_x, shape_y, shape_l;

    forAll(this->nodes_, nodei)
    {
        const vector& ip = nodes_[nodei];
        GaussianQuadrature::calcChebyshev(order_, ip.x(), shape_x);
        GaussianQuadrature::calcChebyshev(order_, ip.y(), shape_y);
        GaussianQuadrature::calcChebyshev
        (
            order_,
            1.0 - ip.x() - ip.y(),
            shape_l
        );

        ni = 0;
        for (label j = 0; j <= order; j++)
        {
            for (label i = 0; i+j <= order; i++)
            {
                T(ni++, nodei) = shape_x[i]*shape_y[j]*shape_l[order-i-j];
            }
        }
    }
    LUscalarMatrix(T).inv(invT_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::triangle::~triangle()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::triangle::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize((order_ + 1)*(order_ + 2)/2, vector::zero);
    G_.setSize(order_*order_);
    E_.setSize(3*P_.size());

    label k = 0;
    for (label j = 0; j <= order_; j++)
    {
        for (label i = 0; i+j <= order_; i++)
        {
            scalar w = cp[i] + cp[j] + cp[order_-i-j];
            P_[k].x() = cp[i]/w;
            P_[k].y() = cp[j]/w;
            k++;
        }
    }

    k = 0;
    label l = 0;
    for (label j = 0; j < order_; j++, k++)
    {
        for (label i = 0; i+j < order_; i++, k++)
        {
            G_[l++] = {k, k+1, k+order_-j+1};
            if (i+j+1 < order_)
            {
                G_[l++] = {k+1, k+order_-j+2, k+order_-j+1};
            }
        }
    }

    label be = 3*order_;
    label ie = 0;
    for (label k = 0; k < order_; k++)
    {
        label& ei = (k == 0 ? be : ie);
        label j = k*(order_ + 1) - (k - 1)*k/2;
        for (label i = 0; i+k < order_; i++)
        {
            E_[ei++] = {j, j+1};
            j++;
        }
    }
    for (label k = order_; k > 0; k--)
    {
        label& ei = (k == order_ ? be : ie);
        label j = k;
        for (label i = 0; i < k; i++)
        {
            E_[ei++] = {j, j+order_-1};
            j += order_-1;
        }
    }
    for (label k = 0; k < order_; k++)
    {
        label& ei = (k == 0 ? be : ie);
        label j = k;
        for (label i = 0; i+k < order_; i++)
        {
            E_[ei++] = {j, j+order_-i+1};
            j += order_-i+1;
        }
    }
}

Foam::tmp<Foam::pointField> Foam::finiteElements::triangle::getNodes
(
    const vector& p0,
    const vector& p1,
    const vector& p2
) const
{
    tmp<pointField> tnodes(new pointField(this->nNodes()));
    pointField& nodes = tnodes.ref();

    vector d10(p1 - p0);
    vector d21(p2 - p1);
    vector d02(p0 - p2);

    const CmptList<vector> x(seg_.nodes(), vector::X);

    nodes[0] = p0;
    nodes[1] = p1;
    nodes[2] = p2;

    label ni = 3;
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = p0 + x[i]*d10;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = p1 + x[i]*d21;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = p2 + x[i]*d02;
    }

    for (label j = 1; j < order_; j++)
    {
        for (label i = 1; i+j < order_; i++)
        {
//             scalar w = 1.0 - (x[i] + x[j]);
            nodes[ni++] = p0 + d10*x[i] - d02*x[j];
        }
    }
    return tnodes;
}


Foam::tmp<Foam::pointField> Foam::finiteElements::triangle::getNodes
(
    const List<vector>& verts
) const
{
    return getNodes(verts[0], verts[1], verts[2]);
}


Foam::scalarList Foam::finiteElements::triangle::calcShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarList(this->nNodes(), 1.0);
    }

    Field<scalar> u(this->nNodes());
    List<scalar> shape_x, shape_y, shape_l;

    GaussianQuadrature::calcChebyshev(order_, pt.x(), shape_x);
    GaussianQuadrature::calcChebyshev(order_, pt.y(), shape_y);
    GaussianQuadrature::calcChebyshev(order_, 1.0 - pt.x() - pt.y(), shape_l);

    label ni = 0;
    forAll(shape_y, j)
    {
        for (label i = 0; i+j <= order_; i++)
        {
            u[ni++] = shape_x[i]*shape_y[j]*shape_l[order_-i-j];
        }
    }
    return invT_*u;
}


Foam::scalarRectangularMatrix Foam::finiteElements::triangle::calcDShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarRectangularMatrix(this->nNodes(), this->nDims(), 0.0);
    }

    scalarRectangularMatrix du(this->nNodes(), this->nDims());
    List<scalar> shape_x, shape_y, shape_l;
    List<scalar> dshape_x, dshape_y, dshape_l;

    GaussianQuadrature::calcChebyshev(order_, pt.x(), shape_x, dshape_x);
    GaussianQuadrature::calcChebyshev(order_, pt.y(), shape_y, dshape_y);
    GaussianQuadrature::calcChebyshev
    (
        order_,
        1.0 - pt.x() - pt.y(),
        shape_l,
        dshape_l
    );

    label ni = 0;
    forAll(shape_y, j)
    {
        for (label i = 0; i+j <= order_; i++)
        {
            const label k = order_-i-j;
            du(ni, 0) =
                (
                    dshape_x[i]*shape_l[k]
                  - shape_x[i]*dshape_l[k]
                )*shape_y[j];
            du(ni, 1) =
                (
                    dshape_y[j]*shape_l[k]
                  - shape_y[j]*dshape_l[k]
                )*shape_x[i];
            ni++;
        }
    }
    return invT_*du;
}


void Foam::finiteElements::triangle::calcDShape
(
    const vector& pt,
    scalarList& shape,
    scalarRectangularMatrix& dshape
) const
{
    shape.setSize(this->nNodes());
    dshape.setSize(this->nNodes(), this->nDims());
    if (order_ == 0)
    {
        shape = 1.0;
        dshape = Zero;
        return;
    }

    scalarRectangularMatrix du(this->nNodes(), this->nDims());
    Field<scalar> u(this->nNodes());
    List<scalar> shape_x, shape_y, shape_l;
    List<scalar> dshape_x, dshape_y, dshape_l;

    GaussianQuadrature::calcChebyshev(order_, pt.x(), shape_x, dshape_x);
    GaussianQuadrature::calcChebyshev(order_, pt.y(), shape_y, dshape_y);
    GaussianQuadrature::calcChebyshev
    (
        order_,
        1.0 - pt.x() - pt.y(),
        shape_l,
        dshape_l
    );

    label ni = 0;
    forAll(shape_y, j)
    {
        for (label i = 0; i+j <= order_; i++)
        {
            const label k = order_-i-j;
            u[ni] = shape_x[i]*shape_y[j]*shape_l[order_-i-j];;
            du(ni, 0) =
                (
                    dshape_x[i]*shape_l[k]
                  - shape_x[i]*dshape_l[k]
                )*shape_y[j];
            du(ni, 1) =
                (
                    dshape_y[j]*shape_l[k]
                  - shape_y[j]*dshape_l[k]
                )*shape_x[i];
            ni++;
        }
    }

    shape = invT_*u;
    dshape = invT_*du;
}


Foam::label Foam::finiteElements::triangle::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 5; // TRIANGLE
        case 2:
            return 22; // QUADRATIC_TRIANGLE
        default:
            return 69; // LAGRANGE_TRIANGLE
    }
}

Foam::label Foam::finiteElements::triangle::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 1;
        case 2:
            return 6;
        case 3:
            return 21;
        case 4:
            return 23;
        case 5:
            return 25;
        default:
            return -1;
    }
}


void Foam::finiteElements::triangle::vtkData
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
