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

#include "hexahedronFiniteElement.H"
#include "CmptList.H"
#include "addToRunTimeSelectionTable.H"
#include "addToRunTimeSelectionMap.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(hexahedron, 0);
    addToRunTimeSelectionTable(finiteElement, hexahedron, type);
    addNamedToRunTimeSelectionTable(finiteElement, hexahedron, type, hex);

    addFE(hex8); addFEMap(hex8, msh, 5);
    addFE(hex27); addFEMap(hex27, msh, 12);
    addFE(hex64); addFEMap(hex64, msh, 92);
    addFE(hex125); addFEMap(hex125, msh, 93);
    addFE(hex216); addFEMap(hex216, msh, 94);
    addFE(hex343); addFEMap(hex343, msh, 95);
    addFE(hex512); addFEMap(hex512, msh, 96);
    addFE(hex729); addFEMap(hex729, msh, 97);
    addFE(hex1000); addFEMap(hex1000, msh, 98);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::hexahedron::hexahedron(const label order)
:
    FiniteElement<ElementType::HEX>
    (
        order,
        (order + 1)*(order + 1)*(order + 1)
    ),
    seg_
    (
        dynamicCast<const segment>
        (
            *finiteElement::getRefFiniteElement(ElementType::SEG, order)
        )
    ),
    map_(finiteElement::dofMap(3, order_))
{
    label ni = 0;
    const CmptList<vector> x(seg_.nodes(), vector::X);

    forAll(x, k)
    {
        forAll(x, j)
        {
            forAll(x, i)
            {
                const label mni = map_[ni++];
                nodes_[mni].x() = x[i];
                nodes_[mni].y() = x[j];
                nodes_[mni].z() = x[k];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::hexahedron::~hexahedron()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::finiteElements::hexahedron::createUniformGeometry() const
{
    const scalarList cp(this->linspace(order_ + 1));
    P_.setSize((order_ + 1)*(order_ + 1)*(order_ + 1), vector::zero);
    G_.setSize(order_*order_*order_);
    E_.setSize(0);

    label l = 0;
    for (label k = 0; k <= order_; k++)
    {
        for (label j = 0; j <= order_; j++)
        {
            for (label i = 0; i <= order_; i++)
            {
                P_[l].x() = cp[i];
                P_[l].y() = cp[j];
                P_[l].z() = cp[k];
                l++;
            }
        }
    }

    l = 0;
    for (label k = 0; k < order_; k++)
    {
        for (label j = 0; j < order_; j++)
        {
            for (label i = 0; i < order_; i++)
            {
                G_[l++] =
                {
                    (i) + ((j) + (k)*(order_+1))*(order_+1),
                    (i+1) + ((j) + (k)*(order_+1))*(order_+1),
                    (i+1) + ((j+1) + (k)*(order_+1))*(order_+1),
                    (i) + ((j+1) + (k)*(order_+1))*(order_+1),
                    (i) + ((j) + (k+1)*(order_+1))*(order_+1),
                    (i+1) + ((j) + (k+1)*(order_+1))*(order_+1),
                    (i+1) + ((j+1) + (k+1)*(order_+1))*(order_+1),
                    (i) + ((j+1) + (k+1)*(order_+1))*(order_+1)
                };
            }
        }
    }
}


Foam::tmp<Foam::pointField> Foam::finiteElements::hexahedron::getNodes
(
    const List<vector>& verts
) const
{
    tmp<pointField> tnodes(new pointField(this->nNodes()));
    pointField& nodes = tnodes.ref();

    const vector& p000(verts[0]);
    const vector& p100(verts[1]);
    const vector& p110(verts[2]);
    const vector& p010(verts[3]);
    const vector& p001(verts[4]);
    const vector& p101(verts[5]);
    const vector& p111(verts[6]);
    const vector& p011(verts[7]);

    const CmptList<vector> x(seg_.nodes(), vector::X);
    label I = 0;
    forAll(x, k)
    {
        const scalar wz = x[k];
        forAll(x, j)
        {
            const scalar wy = x[j];
            forAll(x, i)
            {
                scalar wx = x[i];
                nodes[map_[I++]] =
                    (1.0 - wz)
                   *(
                        (1.0 - wy)*((1.0 - wx)*p000 + wx*p100)
                      + wy*((1.0 - wx)*p010 + wx*p110)
                    )
                  + wz
                   *(
                        (1.0 - wy)*((1.0 - wx)*p001 + wx*p101)
                      + wy*((1.0 - wx)*p011 + wx*p111)
                    );
            }
        }
    }
    return tnodes;
}


Foam::scalarList Foam::finiteElements::hexahedron::calcShape
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
    List<scalar> shape_z(seg_.calcShape(pt.z()));

    label ni = 0;
    forAll(shape_z, k)
    {
        forAll(shape_y, j)
        {
            forAll(shape_x, i)
            {
                shape[map_[ni++]] = shape_x[i]*shape_y[j]*shape_z[k];
            }
        }
    }
    return shape;
}


Foam::scalarRectangularMatrix Foam::finiteElements::hexahedron::calcDShape
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
    List<scalar> shape_z, dshape_z;

    seg_.calcDShape(pt.x(), shape_x, dshape_x);
    seg_.calcDShape(pt.y(), shape_y, dshape_y);
    seg_.calcDShape(pt.z(), shape_z, dshape_z);

    label ni = 0;
    forAll(shape_z, k)
    {
        forAll(shape_y, j)
        {
            forAll(shape_x, i)
            {
                const label mni = map_[ni++];
                dshape(mni, 0) = dshape_x[i]*shape_y[j]*shape_z[k];
                dshape(mni, 1) = shape_x[i]*dshape_y[j]*shape_z[k];
                dshape(mni, 2) = shape_x[i]*shape_y[j]*dshape_z[k];
            }
        }
    }
    return dshape;
}


Foam::label Foam::finiteElements::hexahedron::vtkIndex() const
{
    switch (order_)
    {
        case 1:
            return 12; // CUBE
        case 2:
            return 29; // TRIQUARATIC_CUBE
        default:
            return 72; // LAGRANGE_CUBE
    }
}

Foam::label Foam::finiteElements::hexahedron::mshIndex() const
{
    switch (order_)
    {
        case 1:
            return 5;
        case 2:
            return 12;
        default:
            return -1;
    }
}

void Foam::finiteElements::hexahedron::vtkData
(
    labelList& data,
    label& start,
    const labelList& labels
) const
{
    switch (order_)
    {
        case 0:
        case 1:
        {
            forAll(labels, i)
            {
                data[start++] = labels[i];
            }
            return;
        }
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
            if (!vtkConnectivity_.size())
            {
                vtkConnectivity_ =  vtkElementConnectivity(type_, order_);
            }
            forAll(labels, i)
            {
                data[start++] = labels[vtkConnectivity_[i]];
            }
            return;
        }
    }
}


// ************************************************************************* //
