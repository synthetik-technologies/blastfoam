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

#include "quadrilateralFiniteElement.H"
#include "CmptList.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(quadrilateral, 0);

    addToRunTimeSelectionTable(finiteElement, quadrilateral, type);
    addNamedToRunTimeSelectionTable(finiteElement, quadrilateral, type, quad);

    addFE(quad4); addFEMap(quad4, msh, 3);
    addFE(quad9); addFEMap(quad9, msh, 10);
    addFE(quad16); addFEMap(quad16, msh, 36);
    addFE(quad25); addFEMap(quad25, msh, 37);
    addFE(quad36); addFEMap(quad36, msh, 38);
    addFE(quad49); addFEMap(quad49, msh, 47);
    addFE(quad64); addFEMap(quad64, msh, 48);
    addFE(quad81); addFEMap(quad81, msh, 49);
    addFE(quad100); addFEMap(quad100, msh, 50);
    addFE(quad121); addFEMap(quad121, msh, 51);


//     typedef ShellFiniteElement<quadrilateral> quadrilateralShell;
//     defineTemplateTypeNameAndDebug(quadrilateralShell, 0);
//
//     addToRunTimeSelectionTable
//     (
//         shellFiniteElement,
//         quadrilateralShell,
//         type
//     );
//     addNamedToRunTimeSelectionTable
//     (
//         shellFiniteElement,
//         quadrilateralShell,
//         type,
//         quad
//     );
//     addShellFE(quad4);
//     addShellFE(quad9);
//     addShellFE(quad16);
//     addShellFE(quad25);
//     addShellFE(quad36);
//     addShellFE(quad49);
//     addShellFE(quad64);
//     addShellFE(quad81);
//     addShellFE(quad100);
//     addShellFE(quad121);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::quadrilateral::quadrilateral(const label order)
:
    FiniteElement<ElementType::QUAD>
    (
        order,
        (order + 1)*(order + 1)
    ),
    seg_
    (
        dynamicCast<const segment>
        (
            *finiteElement::getRefFiniteElement(ElementType::SEG, order))
    ),
    map_(finiteElement::dofMap(2, order_))
{
    label ni = 0;
    const CmptList<vector> x(seg_.nodes(), vector::X);

    forAll(x, j)
    {
        forAll(x, i)
        {
            const label mni = map_[ni++];
            nodes_[mni].x() = x[i];
            nodes_[mni].y() = x[j];
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::quadrilateral::~quadrilateral()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::quadrilateral::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize((order_ + 1)*(order_ + 1), vector::zero);
    G_.setSize(order_*order_);
    E_.setSize(2*order_*(order_ + 1));

    label k = 0;
    for (label j = 0; j <= order_; j++)
    {
        for (label i = 0; i <= order_; i++)
        {
            P_[k].x() = cp[i];
            P_[k].y() = cp[j];
            k++;
        }
    }

    k = 0;
    label l = 0;
    for (label j = 0; j < order_; j++, k++)
    {
        for (label i = 0; i < order_; i++, k++)
        {
            G_[l++] =
            {
                k,
                k+1,
                k+order_+2,
                k+order_+1
            };
        }
    }

    label be = 4*order_;
    label ie = 0;
    for (label k = 0; k <= order_; k++)
    {
        label& ei = (k == 0 || k == order_ ? be : ie);
        for (label i = 0, j = k*(order_+1); i < order_; i++)
        {
            E_[ei++] = {j, j+1};
            j++;
        }
    }
    for (label k = order_; k >= 0; k--)
    {
        label& ei = (k = 0 || k == order_ ? be : ie);
        for (label i = 0, j = k; i < order_; i++, j+= order_+1)
        {
            E_[ei++] = {j, j+order_+1};
        }
    }
}


Foam::tmp<Foam::pointField> Foam::finiteElements::quadrilateral::getNodes
(
    const List<vector>& verts
) const
{
    tmp<pointField> tnodes(new pointField(this->nNodes()));
    pointField& nodes = tnodes.ref();

    const vector& p00(verts[0]);
    const vector& p10(verts[1]);
    const vector& p11(verts[2]);
    const vector& p01(verts[3]);

    const CmptList<vector> x(seg_.nodes(), vector::X);
    label I = 0;
    for (label j = 0; j < x.size(); j++)
    {
        scalar wy = x[j];
        for (label i = 0; i < x.size(); i++)
        {
            scalar wx = x[i];
            nodes[map_[I++]] =
                (1.0 - wy)*((1.0 - wx)*p00 + wx*p10)
              + wy*((1.0 - wx)*p01 + wx*p11);
        }
    }
    return tnodes;
}


Foam::scalarList Foam::finiteElements::quadrilateral::calcShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarList(this->nNodes(), 1.0);
    }

    scalarList shape(this->nNodes());
    List<scalar> shape_x(seg_.calcShape(pt.x()));
    List<scalar> shape_y(seg_.calcShape(pt.y()));

    label ni = 0;
    forAll(shape_y, j)
    {
        forAll(shape_x, i)
        {
            shape[map_[ni++]] = shape_x[i]*shape_y[j];
        }
    }
    return shape;
}


Foam::scalarRectangularMatrix Foam::finiteElements::quadrilateral::calcDShape
(
    const vector& pt
) const
{
    if (order_ == 0)
    {
        return scalarRectangularMatrix(this->nNodes(), this->nDims(), 0.0);
    }

    scalarRectangularMatrix dshape(this->nNodes(), this->nDims());
    List<scalar> shape_x, dshape_x;
    List<scalar> shape_y, dshape_y;

    seg_.calcDShape(pt.x(), shape_x, dshape_x);
    seg_.calcDShape(pt.y(), shape_y, dshape_y);

    label ni = 0;
    forAll(shape_y, j)
    {
        forAll(shape_x, i)
        {
            const label mni = map_[ni++];
            dshape(mni, 0) = dshape_x[i]*shape_y[j];
            dshape(mni, 1) = shape_x[i]*dshape_y[j];
        }
    }
    return dshape;
}


Foam::label Foam::finiteElements::quadrilateral::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 9; // SQUARE
        case 2:
            return 28; // BIQUADRATIC_SQUARE
        default:
            return 70; // LAGRANGE_SQUARE
    }
}

Foam::label Foam::finiteElements::quadrilateral::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 3;
        case 2:
            return 10;
        default:
            return -1;
    }
}


void Foam::finiteElements::quadrilateral::vtkData
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
