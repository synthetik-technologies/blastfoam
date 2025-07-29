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

#include "prismFiniteElement.H"
#include "CmptList.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(prism, 0);

    addToRunTimeSelectionTable(finiteElement, prism, type);

    addFE(wedge6); addFEMap(wedge6, msh, 6);
    addFE(wedge18); addFEMap(wedge18, msh, 13);
    addFE(wedge40); addFEMap(wedge40, msh, 90);
    addFE(wedge75); addFEMap(wedge75, msh, 91);
    addFE(wedge126); addFEMap(wedge126, msh, 106);
    addFE(wedge196); addFEMap(wedge196, msh, 107);
    addFE(wedge288); addFEMap(wedge288, msh, 108);
    addFE(wedge405); addFEMap(wedge405, msh, 109);
    addFE(wedge550); addFEMap(wedge550, msh, 110);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::prism::prism(const label order)
:
    FiniteElementBase<ElementType::PRISM>
    (
        order,
        (order + 1)*(order + 1)*(order + 2)/2.0
    ),
    seg_
    (
        dynamicCast<const segment>
        (
            *finiteElement::getRefFiniteElement(ElementType::SEG, order))
    ),
    tri_
    (
        dynamicCast<const triangle>
        (
            *finiteElement::getRefFiniteElement(ElementType::TRI, order))
    )
{
    const CmptList<vector> x(seg_.nodes(), vector::X);
    const List<vector>& triNodes = tri_.nodes();

    label ni = 0;
    forAll(x, i)
    {
        forAll(triNodes, j)
        {
            nodes_[ni][0] = triNodes[j][0];
            nodes_[ni][1] = triNodes[j][1];
            nodes_[ni][2] = x[i];
            ni++;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::prism::~prism()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::prism::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize((order_ + 1)*(order_ + 1)*(order_ + 2)/2, vector::zero);
    G_.setSize(order_*order_*order_);
    E_.setSize(0);

    for (label k = 0, m = 0; k <= order_; k++)
    {
        for (label j = 0; j <= order_; j++)
        {
            for (label i = 0; i+j <= order_; i++)
            {
                scalar w = cp[i] + cp[j] + cp[order_-i-j];
                P_[m].x() = cp[i]/w;
                P_[m].y() = cp[j]/w;
                P_[m].z() = cp[k];
                m++;
            }
        }
    }

    for (label k = 0; k < order_; k++)
    {
        for (label j = 0, l = 0; j < order_; j++)
        {
            for (label i = 0; i < order_; i++)
            {
                G_[l++] =
                {
                    l + (k)*(order_+1)*(order_+2)/2,
                    l + 1 + (k)*(order_+1)*(order_+2)/2,
                    l - j + (2 + (k)*(order_+2))*(order_+1)/2,
                    l + (k+1)*(order_+1)*(order_+2)/2,
                    l + 1 + (k+1)*(order_+1)*(order_+2)/2,
                    l - j + (2 + (k+1)*(order_+2))*(order_+1)/2
                };
                if (i+j+1 < order_)
                {
                    G_[l++] =
                    {
                        l + 1 + (k)*(order_+1)*(order_+2)/2,
                        l - j + (2 + (k)*(order_+1))*(order_+2)/2,
                        l - j + (2 + (k)*(order_+2))*(order_+1)/2,
                        l + 1 + (k+1)*(order_+1)*(order_+2)/2,
                        l - j + (2 + (k+1)*(order_+1))*(order_+2)/2,
                        l - j + (2 + (k+1)*(order_+2))*(order_+1)/2,
                    };
                }
            }
        }
    }
}


Foam::tmp<Foam::pointField> Foam::finiteElements::prism::getNodes
(
    const List<vector>& verts
) const
{
    tmp<pointField> tnodes(new pointField(this->nNodes()));
    pointField& nodes = tnodes.ref();

    const vector& p00(verts[0]);
    const vector& p10(verts[1]);
    const vector& p20(verts[2]);

    const vector& p01(verts[3]);
    const vector& p11(verts[4]);
    const vector& p21(verts[5]);

    const CmptList<vector> x(seg_.nodes(), vector::X);
    Field<vector> t1(tri_.getNodes(p00, p10, p20));
    Field<vector> t2(tri_.getNodes(p01, p11, p21));

    label ni = 0;
    forAll(x, i)
    {
        forAll(t1, j)
        {
            nodes[ni++] = (1.0 - x[i])*t1[j] + x[i]*t2[j];
        }
    }
    return tnodes;
}


Foam::scalarList Foam::finiteElements::prism::calcShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarList(this->nNodes(), 1.0);
    }

    scalarList shape(this->nNodes());
    List<scalar> shape_t(tri_.calcShape(pt));
    List<scalar> shape_s(seg_.calcShape(pt.z()));

    label ni = 0;
    forAll(shape_s, i)
    {
        forAll(shape_t, j)
        {
            shape[ni++] = shape_s[i]*shape_t[j];
        }
    }
    return shape;
}


Foam::scalarRectangularMatrix Foam::finiteElements::prism::calcDShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarRectangularMatrix(this->nNodes(), this->nDims(), 0.0);
    }

    scalarRectangularMatrix dshape(this->nNodes(), this->nDims());
    List<scalar> shape_t;
    scalarRectangularMatrix dshape_t;
    List<scalar> shape_s, dshape_s;

    tri_.calcDShape(pt, shape_t, dshape_t);
    seg_.calcDShape(pt.z(), shape_s, dshape_s);

    label ni = 0;
    forAll(shape_s, i)
    {
        forAll(shape_t, j)
        {
            dshape(ni, 0) = dshape_t(j, 0)*shape_s[i];
            dshape(ni, 1) = dshape_t(j, 1)*shape_s[i];
            dshape(ni, 2) = shape_t[j]*dshape_s[i];
            ni++;
        }
    }
    return dshape;
}


Foam::label Foam::finiteElements::prism::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 13; // PRISM
        case 2:
            return 32; // BI_QUADRTIC_PRISM
        default:
            return 73; // LAGRANGE_PRISM
    }
}

Foam::label Foam::finiteElements::prism::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 6;
        case 2:
            return 13;
        default:
            return -1;
    }
}


void Foam::finiteElements::prism::vtkData
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
        case 2:
        {
            data[start++] = labels[0];
            data[start++] = labels[2];
            data[start++] = labels[1];
            data[start++] = labels[3];
            data[start++] = labels[5];
            data[start++] = labels[4];
            data[start++] = labels[8];
            data[start++] = labels[7];
            data[start++] = labels[6];
            data[start++] = labels[11];
            data[start++] = labels[10];
            data[start++] = labels[9];
            data[start++] = labels[12];
            data[start++] = labels[14];
            data[start++] = labels[13];
            data[start++] = labels[17];
            data[start++] = labels[16];
            data[start++] = labels[15];
            return;
        }
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
