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

#include "pyramidFiniteElement.H"
#include "CmptList.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(pyramid, 0);

    addToRunTimeSelectionTable(finiteElement, pyramid, type);
    addNamedToRunTimeSelectionTable(finiteElement, pyramid, type, pyr);

    addFE(pyramid5); addFEMap(pyramid5, msh, 7);
//     addFE(pyramid14); addFEMap(pyramid14, msh, 14);
//     addFE(pyramid30); addFEMap(pyramid30, msh, 118);
//     addFE(pyramid55); addFEMap(pyramid55, msh, 119);
//     addFE(pyramid91); addFEMap(pyramid91, msh, 120);
//     addFE(pyramid140); addFEMap(pyramid140, msh, 121);
//     addFE(pyramid204); addFEMap(pyramid204, msh, 122);
//     addFE(pyramid285); addFEMap(pyramid285, msh, 123);
//     addFE(pyramid385); addFEMap(pyramid385, msh, 124);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::pyramid::pyramid(const label order)
:
    FiniteElement<ElementType::PYR>
    (
        order,
        5
    )
{
    if (order > 1)
    {
        FatalErrorInFunction
            << "Only 1st order pyramids are supported" << endl
            << abort(FatalError);
    }
    nodes_[0].x() = 0.0;
    nodes_[0].y() = 0.0;
    nodes_[0].z() = 0.0;

    nodes_[1].x() = 1.0;
    nodes_[1].y() = 0.0;
    nodes_[1].z() = 0.0;

    nodes_[2].x() = 1.0;
    nodes_[2].y() = 1.0;
    nodes_[2].z() = 0.0;

    nodes_[3].x() = 0.0;
    nodes_[3].y() = 1.0;
    nodes_[3].z() = 0.0;

    nodes_[4].x() = 0.0;
    nodes_[4].y() = 0.0;
    nodes_[4].z() = 1.0;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::pyramid::~pyramid()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::pyramid::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize((order_ + 1)*(order_ + 1)*(2*order_ + 3)/6, vector::zero);
    G_.setSize(order_*(2*order_ - 1)*(2*order_ + 1)/3);
    E_.setSize(0);

    for (label k = 0, m = 0; k <= order_; k++)
    {
        const scalarList cpk(this->linspace(order_ - k));
        for (label j = 0; j+k <= order_; j++)
        {
            for (label i = 0; i+k <= order_; i++)
            {
//                 P_[m].x() = (order_ > k ? scalar(i)/(scalar(order_ - k) : 0.0);
//                 P_[m].y() = (order_ > k ? scalar(j)/(scalar(order_ - k) : 0.0);
//                 P_[m].z() = scalar(k)/scalar(order_);

                P_[m].x() = cpk[i]*(1.0 - cp[k]);
                P_[m].y() = cpk[j]*(1.0 - cp[k]);
                P_[m].z() = cp[k];
                m++;
            }
        }
    }

    for (label k = 0, m = 0; k < order_; k++)
    {
        label lk = k*(k*(2*k - 6*order_ - 9) + 6*order_*(order_ + 3) + 13)/6;
        label lkp1 =
            (k + 1)*(k*(2*k - 6*order_ - 5) + 6*order_*(order_ + 2) + 6)/6;
        for (label j = 0; j+k < order_; j++)
        {
            for (label i = 0; i+k < order_; i++)
            {
                G_[m++] =
                {
                    lk + j*(order_ - k + 1) + i,
                    lk + j*(order_ - k + 1) + i + 1,
                    lk + (j + 1)*(order_ - k + 1) + i + 1,
                    lk + (j + 1)*(order_ - k + 1) + i,
                    lkp1 + j*(order_ - k) + i,
                };

            }
        }
        for (label j = 0; j+k+1 < order_; j++)
        {
            for (label i = 0; i+k+1 < order_; i++)
            {
                G_[m++] =
                {
                    lkp1 + j*(order_ - k) + i,
                    lkp1 + (j + 1)*(order_ - k) + i,
                    lkp1 + (j + 1)*(order_ - k) + i + 1,
                    lkp1 + j*(order_ - k) + i + 1,
                    lk + (j + 1)*(order_ - k + 1) + i + 1,
                };

            }
        }
        for (label j = 0; j+k < order_; j++)
        {
            for (label i = 0; i+k+1 < order_; i++)
            {
                G_[m++] =
                {
                    lk + j*(order_ - k + 1) + i + 1,
                    lk + (j + 1)*(order_ - k + 1) + i + 1,
                    lkp1 + j*(order_ - k) + i,
                    lkp1 + j*(order_ - k) + i + 1
                };

            }
        }
        for (label j = 0; j+k+1 < order_; j++)
        {
            for (label i = 0; i+k < order_; i++)
            {
                G_[m++] =
                {
                    lk + (j + 1)*(order_ - k + 1) + i,
                    lk + (j + 1)*(order_ - k + 1) + i + 1,
                    lkp1 + (j + 1)*(order_ - k) + i,
                    lkp1 + j*(order_ - k) + i
                };

            }
        }
    }
}


Foam::tmp<Foam::pointField> Foam::finiteElements::pyramid::getNodes
(
    const List<vector>& verts
) const
{
    return tmp<pointField>(new pointField(verts));
}


Foam::scalarList Foam::finiteElements::pyramid::calcShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarList(this->nNodes(), 1.0);
    }
Info<<this->nNodes()<<endl;
    scalarList shape(this->nNodes());
    const scalar z = pt.z();
    const scalar oz = 1.0 - z;
    if (oz < 1e-6)
    {
        shape = 0.0;
        shape[4] = 1.0;
    }
    else
    {
        const scalar x = pt.x();
        const scalar y = pt.y();
        const scalar ox = 1.0 - x - z;
        const scalar oy = 1.0 - y - z;
        shape[0] = ox*oy/oz;
        shape[1] = x*oy/oz;
        shape[2] = x*y/oz;
        shape[3] = ox*y/ox;
        shape[4] = z;
    }
    return shape;
}


Foam::scalarRectangularMatrix Foam::finiteElements::pyramid::calcDShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarRectangularMatrix(this->nNodes(), this->nDims(), 0.0);
    }

    scalarRectangularMatrix dshape(this->nNodes(), this->nDims());

    const scalar z = pt.z();
    const scalar oz = 1.0 - z;
    if (oz < 1e-6)
    {
        dshape(0, 0) = -0.5;
        dshape(0, 1) = -0.5;
        dshape(0, 2) = -0.75;

        dshape(1, 0) = 0.5;
        dshape(1, 1) = -0.5;
        dshape(1, 2) = -0.25;

        dshape(2, 0) = 0.5;
        dshape(2, 1) = 0.5;
        dshape(2, 2) = 0.25;

        dshape(3, 0) = -0.5;
        dshape(3, 1) = 0.5;
        dshape(3, 2) = 0.25;

        dshape(4, 0) = 0.0;
        dshape(4, 1) = 0.0;
        dshape(4, 2) = 1.0;
    }
    else
    {
        const scalar x = pt.x();
        const scalar y = pt.y();
        const scalar ox = 1.0 - x - z;
        const scalar oy = 1.0 - y - z;

        dshape(0, 0) = -oy/oz;
        dshape(0, 1) = -ox/oz;
        dshape(0, 2) = x*y/sqr(oz) - 1.0;

        dshape(1, 0) = oy/oz;
        dshape(1, 1) = -x/oz;
        dshape(1, 2) = -x*y/sqr(oz);

        dshape(2, 0) = y/oz;
        dshape(2, 1) = x/oz;
        dshape(2, 2) = x*y/sqr(oz);

        dshape(3, 0) = -y/oz;
        dshape(3, 1) = ox/oz;
        dshape(3, 2) = -x*y/sqr(oz);

        dshape(4, 0) = 0.0;
        dshape(4, 1) = 0.0;
        dshape(4, 2) = 1.0;
    }
    return dshape;
}


Foam::label Foam::finiteElements::pyramid::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 14; // PYRAMID
        case 2:
            return 27; // QUADRATIC_PYRAMID
        default:
            return 74; // LAGRANGE_PYRAMID
    }
}

Foam::label Foam::finiteElements::pyramid::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 7;
        default:
            return -1;
    }
}


void Foam::finiteElements::pyramid::vtkData
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
