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

#include "tetrahedronFiniteElement.H"
#include "triangleFiniteElement.H"
#include "CmptList.H"
#include "GaussianQuadrature.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(tetrahedron, 0);

    addToRunTimeSelectionTable(finiteElement, tetrahedron, type);
    addNamedToRunTimeSelectionTable(finiteElement, tetrahedron, type, tet);

    addFE(tet4); addFEMap(tet4, msh, 4);
    addFE(tet10); addFEMap(tet10, msh, 11);
    addFE(tet20); addFEMap(tet20, msh, 29);
    addFE(tet35); addFEMap(tet35, msh, 30);
    addFE(tet56); addFEMap(tet56, msh, 31);
    addFE(tet84); addFEMap(tet84, msh, 71);
    addFE(tet120); addFEMap(tet120, msh, 72);
    addFE(tet165); addFEMap(tet165, msh, 73);
    addFE(tet220); addFEMap(tet220, msh, 74);
    addFE(tet286); addFEMap(tet286, msh, 75);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::tetrahedron::tetrahedron(const label order)
:
    FiniteElement<ElementType::TET>
    (
        order,
        (order + 1)*(order + 2)*(order + 3)/6
    ),
    seg_
    (
        dynamicCast<const segment>
        (
            *finiteElement::getRefFiniteElement(ElementType::SEG, order))
    ),
    nFaceNodes_(triangle::calcNFaceNodes(order_)),
    invT_()
{
    label ni = 0;

    const CmptList<vector> x(seg_.nodes(), vector::X);
    nodes_[0][0] = x[0];
    nodes_[0][1] = x[0];
    nodes_[0][2] = x[0];

    nodes_[1][0] = x[order];
    nodes_[1][1] = x[0];
    nodes_[1][2] = x[0];

    nodes_[2][0] = x[0];
    nodes_[2][1] = x[order];
    nodes_[2][2] = x[0];

    nodes_[3][0] = x[0];
    nodes_[3][1] = x[0];
    nodes_[3][2] = x[order];

    ni = 4;
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[i];
        nodes_[ni][1] = x[0];
        nodes_[ni][2] = x[0];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[0];
        nodes_[ni][1] = x[i];
        nodes_[ni][2] = x[0];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[0];
        nodes_[ni][1] = x[0];
        nodes_[ni][2] = x[i];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[order-i];
        nodes_[ni][1] = x[i];
        nodes_[ni][2] = x[0];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[order-i];
        nodes_[ni][1] = x[0];
        nodes_[ni][2] = x[i];
        ni++;
    }
    for (label i = 1; i < order; i++)
    {
        nodes_[ni][0] = x[0];
        nodes_[ni][1] = x[order-i];
        nodes_[ni][2] = x[i];
        ni++;
    }

    for (label j = 1; j < order; j++)
    {
        for (label i = 1; i+j < order; i++)
        {
            scalar w = x[i] + x[j] + x[order-i-j];
            nodes_[ni][0] = x[order-i-j]/w;
            nodes_[ni][1] = x[i]/w;
            nodes_[ni][2] = x[j]/w;
            ni++;
        }
    }
    for (label j = 1; j < order; j++)
    {
        for (label i = 1; i+j < order; i++)
        {
            scalar w = x[i] + x[j] + x[order-i-j];
            nodes_[ni][0] = x[0]/w;
            nodes_[ni][1] = x[j]/w;
            nodes_[ni][2] = x[i]/w;
            ni++;
        }
    }
    for (label j = 1; j < order; j++)
    {
        for (label i = 1; i+j < order; i++)
        {
            scalar w = x[i] + x[j] + x[order-i-j];
            nodes_[ni][0] = x[i]/w;
            nodes_[ni][1] = x[0]/w;
            nodes_[ni][2] = x[j]/w;
            ni++;
        }
    }
    for (label j = 1; j < order; j++)
    {
        for (label i = 1; i+j < order; i++)
        {
            scalar w = x[i] + x[j] + x[order-i-j];
            nodes_[ni][0] = x[j]/w;
            nodes_[ni][1] = x[i]/w;
            nodes_[ni][2] = x[0]/w;
            ni++;
        }
    }
    for (label k = 1; k < order; k++)
    {
        for (label j = 1; j+k < order; j++)
        {
            for (label i = 1; i+j+k < order; i++)
            {
                scalar w =
                    x[i] + x[j] + x[k] + x[order-i-j-k];
                nodes_[ni][0] = x[i]/w;
                nodes_[ni][1] = x[j]/w;
                nodes_[ni][2] = x[k]/w;
                ni++;
            }
        }
    }


    invT_.setSize(this->nNodes());

    List<scalar> shape_x(this->nNodes());
    List<scalar> shape_y(this->nNodes());
    List<scalar> shape_z(this->nNodes());
    List<scalar> shape_l(this->nNodes());

    forAll(this->nodes_, nodei)
    {
        const vector& ip = nodes_[nodei];
        GaussianQuadrature::calcChebyshev(order_, ip.x(), shape_x);
        GaussianQuadrature::calcChebyshev(order_, ip.y(), shape_y);
        GaussianQuadrature::calcChebyshev(order_, ip.z(), shape_z);
        GaussianQuadrature::calcChebyshev
        (
            order_,
            1.0 - ip.x() - ip.y() - ip.z(),
            shape_l
        );

        ni = 0;
        for (label k = 0; k <= order; k++)
        {
            for (label j = 0; j+k <= order; j++)
            {
                for (label i = 0; i+j+k <= order; i++)
                {
                    invT_(ni++, nodei) =
                        shape_x[i]
                       *shape_y[j]
                       *shape_z[k]
                       *shape_l[order-i-j-k];
                }
            }
        }
    }

    const label m = invT_.m();
    pivotIndices_.setSize(m);
    LUDecompose(invT_, pivotIndices_);
//     for (label i = 0; i < m; i++)
//     {
//         // Pivot
//         {
//             label piv = i;
//             scalar a = mag(invT_(piv, i));
//             for (label j = i+1; j < m; j++)
//             {
//                 scalar b = mag(invT_(j, i));
//                 if (b > a)
//                 {
//                     a = b;
//                     piv = j;
//                 }
//             }
//             pivotIndices_[i] = piv;
//             if (piv != i)
//             {
//                 for (label j = 0; j < m; j++)
//                 {
//                     Swap(invT_(i, j), invT_(piv, j));
//                 }
//             }
//         }
//
//         if (mag(invT_(i, i)) < small)
//         {
//             FatalErrorInFunction<<abort(FatalError);
//         }
//
//         const scalar raii = 1.0/invT_(i, i);
//         for (label j = i+1; j < m; j++)
//         {
//             invT_(j, i) *= raii;
//         }
//         for (label k = i+1; k < m; k++)
//         {
//             const scalar aik = invT_(i, k);
//             for (label j = i+1; j < m; j++)
//             {
//                 invT_(j, k) -= aik*invT_(j, i);
//             }
//         }
//     }

//     invT_.decompose();
//     Info<<T<<endl;
//     LUscalarMatrix(T).inv(invT_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::tetrahedron::~tetrahedron()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::tetrahedron::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize((order_ + 1)*(order_ + 2)*(order_+3)/6, vector::zero);
    G_.setSize(order_*order_*order_);
    E_.setSize(0);
    labelList v(cp.size()*cp.size()*cp.size());
    label m = 0;

    for (label k = 0; k <= order_; k++)
    {
        for (label j = 0; j+k <= order_; j++)
        {
            for (label i = 0; i+j+k <= order_; i++)
            {
                scalar w = cp[i] + cp[j] + cp[k] + cp[order_-i-j-k];
                P_[m].x() = cp[i]/w;
                P_[m].y() = cp[j]/w;
                P_[m].y() = cp[k]/w;

                v[j + (j+k + (i+j+k)*(order_+1))*(order_+1)] = m;
                m++;
            }
        }
    }
    if (m != P_.size())
    {
        FatalErrorInFunction
            << "Bad points size" << endl
            << abort(FatalError);
    }

    m = 0;
    for (label k = 0; k < order_; k++)
    {
        for (label j = 0; j <= k; j++)
        {
            for (label i = 0; i <= j; i++)
            {
                G_[m++] =
                {
                    v[(i) + ((j) + (k)*(order_+1))*(order_+1)],
                    v[(i) + ((j) + (k+1)*(order_+1))*(order_+1)],
                    v[(i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1)],
                    v[(i) + ((j+1) + (k+1)*(order_+1))*(order_+1)]
                };
                if (j < k)
                {
                    G_[m++] =
                    {
                        v[(i) + ((j) + (k)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1)],
                        v[(i) + ((j+1) + (k)*(order_+1))*(order_+1)],
                        v[(i) + ((j+1) + (k+1)*(order_+1))*(order_+1)]
                    };

                    G_[m++] =
                    {
                        v[(i) + ((j) + (k)*(order_+1))*(order_+1)],
                        v[(i) + ((j+1) + (k)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j+1) + (k)*(order_+1))*(order_+1)]
                    };
                }
                if (i < j)
                {
                    G_[m++] =
                    {
                        v[(i) + ((j) + (k)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j) + (k)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j) + (k+1)*(order_+1))*(order_+1)]
                    };

                    if (j < k)
                    {
                        G_[m++] =
                        {
                            v[(i) + ((j) + (k)*(order_+1))*(order_+1)],
                            v[(i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1)],
                            v[(i+1) + ((j) + (k)*(order_+1))*(order_+1)],
                            v[(i+1) + ((j+1) + (k)*(order_+1))*(order_+1)]
                        };
                    }
                    G_[m++] =
                    {
                        v[(i) + ((j) + (k)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1)],
                        v[(i) + ((j) + (k+1)*(order_+1))*(order_+1)],
                        v[(i+1) + ((j) + (k+1)*(order_+1))*(order_+1)]
                    };
                }
            }
        }
    }
    if (m != G_.size())
    {
        FatalErrorInFunction
            << "Bad geometry size" << endl
            << abort(FatalError);
    }
}


Foam::tmp<Foam::pointField> Foam::finiteElements::tetrahedron::getNodes
(
    const List<vector>& verts

) const
{
    tmp<pointField> tnodes(new pointField(this->nNodes()));
    pointField& nodes = tnodes.ref();

    const vector& p0(verts[0]);
    const vector& p1(verts[1]);
    const vector& p2(verts[2]);
    const vector& p3(verts[3]);

    const CmptList<vector> x(seg_.nodes(), vector::X);

    label ni = 4;
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = (1.0 - x[i])*p0 + x[i]*p1;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = (1.0 - x[i])*p0 + x[i]*p2;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = (1.0 - x[i])*p0 + x[i]*p3;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = (1.0 - x[i])*p1 + x[i]*p2;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = (1.0 - x[i])*p1 + x[i]*p3;
    }
    for (label i = 1; i < order_; i++)
    {
        nodes[ni++] = (1.0 - x[i])*p2 + x[i]*p3;
    }

    for (label j = 1; j < order_; j++)
    {
        for (label i = 1; i+j < order_; i++)
        {
            scalar w = x[i] + x[j] + x[order_-i-j];
            nodes[ni++] =
                p0*x[order_-i-j]/w
              + p1*x[i]/w
              + p2*x[j]/w;
        }
    }
    for (label j = 1; j < order_; j++)
    {
        for (label i = 1; i+j < order_; i++)
        {
            scalar w = x[i] + x[j] + x[order_-i-j];
            nodes[ni++] =
                p0*x[order_-i-j]/w
              + p1*x[i]/w
              + p3*x[j]/w;
        }
    }
    for (label j = 1; j < order_; j++)
    {
        for (label i = 1; i+j < order_; i++)
        {
            scalar w = x[i] + x[j] + x[order_-i-j];
            nodes[ni++] =
                p0*x[order_-i-j]/w
              + p2*x[i]/w
              + p3*x[j]/w;
        }
    }
    for (label j = 1; j < order_; j++)
    {
        for (label i = 1; i+j < order_; i++)
        {
            scalar w = x[i] + x[j] + x[order_-i-j];
            nodes[ni++] =
                p1*x[order_-i-j]/w
              + p2*x[i]/w
              + p3*x[j]/w;
        }
    }
    for (label k = 1; k < order_; k++)
    {
        for (label j = 1; j+k < order_; j++)
        {
            for (label i = 1; i+j+k < order_; i++)
            {
                nodes[ni++] =
                    (1.0 - x[i] - x[j] - x[k])*p0
                  + p1*x[i]
                  + p2*x[j]
                  + p3*x[k];
            }
        }
    }
    return tnodes;
}


Foam::scalarList Foam::finiteElements::tetrahedron::calcShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarList(this->nNodes(), 1.0);
    }

    Field<scalar> u(this->nNodes());
    List<scalar> shape_x, shape_y, shape_z, shape_l;
    GaussianQuadrature::calcChebyshev(order_, pt.x(), shape_x);
    GaussianQuadrature::calcChebyshev(order_, pt.y(), shape_y);
    GaussianQuadrature::calcChebyshev(order_, pt.z(), shape_z);
    GaussianQuadrature::calcChebyshev
    (
        order_,
        1.0 - pt.x() - pt.y() - pt.z(),
        shape_l
    );

    label ni = 0;
    for (label k = 0; k <= order_; k++)
    {
        for (label j = 0; j+k <= order_; j++)
        {
            for (label i = 0; i+j+k <= order_; i++)
            {
                u[ni++] =
                    shape_x[i]
                   *shape_y[j]
                   *shape_z[k]
                   *shape_l[order_-i-j-k];
            }
        }
    }

    LUBacksubstitute(invT_, pivotIndices_, u);
    return u;
}


Foam::scalarRectangularMatrix Foam::finiteElements::tetrahedron::calcDShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarRectangularMatrix(this->nNodes(), this->nDims(), 0.0);
    }

    List<scalarField> du(3, scalarField(this->nNodes()));
    List<scalar> shape_x, shape_y, shape_z, shape_l;
    List<scalar> dshape_x, dshape_y, dshape_z, dshape_l;

    GaussianQuadrature::calcChebyshev(order_, pt.x(), shape_x, dshape_x);
    GaussianQuadrature::calcChebyshev(order_, pt.y(), shape_y, dshape_y);
    GaussianQuadrature::calcChebyshev(order_, pt.z(), shape_z, dshape_z);
    GaussianQuadrature::calcChebyshev
    (
        order_,
        1.0 - pt.x() - pt.y() - pt.z(),
        shape_l,
        dshape_l
    );

    label ni = 0;
    for (label k = 0; k <= order_; k++)
    {
        for (label j = 0; j+k <= order_; j++)
        {
            for (label i = 0; i+j+k <= order_; i++)
            {
                label l = order_-i-j-k;
                du[0][ni] =
                    (
                        dshape_x[i]*shape_l[l]
                      - shape_x[i]*dshape_l[l]
                    )*shape_y[j]*shape_z[k];
                du[1][ni] =
                    (
                        dshape_y[j]*shape_l[l]
                      - shape_y[j]*dshape_l[l]
                    )*shape_x[i]*shape_z[k];
                du[2][ni] =
                    (
                        dshape_z[k]*shape_l[l]
                      - shape_z[k]*dshape_l[l]
                    )*shape_x[i]*shape_y[j];
                ni++;
            }
        }
    }

    LUBacksubstitute(invT_, pivotIndices_, du[0]);
    LUBacksubstitute(invT_, pivotIndices_, du[1]);
    LUBacksubstitute(invT_, pivotIndices_, du[2]);
    scalarRectangularMatrix dshape(this->nNodes(), this->nDims());
    for (label i = 0; i < invT_.m(); i++)
    {
        dshape(i, 0) = du[0][i];
        dshape(i, 1) = du[1][i];
        dshape(i, 2) = du[2][i];
    }

    return dshape;
}


Foam::label Foam::finiteElements::tetrahedron::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 10; // TETRAHEDRON
        case 2:
            return 24; // QUADRATIC_TETRAHEDRON
        default:
            return 71; // LAGRANGE_TETRAHEDRON
    }
}

Foam::label Foam::finiteElements::tetrahedron::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 4;
        case 2:
            return 11;
        case 3:
            return 29;
        case 4:
            return 30;
        case 5:
            return 31;
        default:
            return -1;
    }
}


void Foam::finiteElements::tetrahedron::vtkData
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
            data[start++] = labels[1];
            data[start++] = labels[2];
            data[start++] = labels[3];
            data[start++] = labels[4];
            data[start++] = labels[5];
            data[start++] = labels[6];
            data[start++] = labels[7];
            data[start++] = labels[8];
            data[start++] = labels[9];
            data[start++] = labels[10];
            data[start++] = labels[11];
            data[start++] = labels[12];
            data[start++] = labels[13];
            data[start++] = labels[14];
            data[start++] = labels[15];
            data[start++] = labels[16];
            data[start++] = labels[17];
            data[start++] = labels[18];
            data[start++] = labels[19];
            data[start++] = labels[24];
            data[start++] = labels[22];
            data[start++] = labels[21];
            data[start++] = labels[23];
            data[start++] = labels[20];
            data[start++] = labels[25];
            data[start++] = labels[26];
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
