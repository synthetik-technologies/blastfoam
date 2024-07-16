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

#include "quadrilateralShellFiniteElements.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(quadrilateral4Shell1, 0);
    defineTypeNameAndDebug(quadrilateral4Shell2, 0);
    defineTypeNameAndDebug(quadrilateral4Shell3, 0);

    addToRunTimeSelectionTable(shellFiniteElement, quadrilateral4Shell1, type);
    addToRunTimeSelectionTable(shellFiniteElement, quadrilateral4Shell2, type);
    addToRunTimeSelectionTable(shellFiniteElement, quadrilateral4Shell3, type);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::finiteElements::quadrilateral4Shell::quadrilateral4Shell()
:
    FiniteElement<ElementType::QUAD>(1, 4)
{
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
}


Foam::finiteElements::quadrilateral4Shell1::quadrilateral4Shell1()
{
    ir_.setSize(4);

    const scalar a = 1.0/sqrt(3.0);
    const scalar na = (1.0 - a)*0.5;
    const scalar pa = 1.0 - na;

    ir_[0].x() = na;
    ir_[0].y() = na;
    ir_[0].z() = 0.0;
    ir_[0].w() = 0.25;

    ir_[1].x() = pa;
    ir_[1].y() = na;
    ir_[1].z() = 0.0;
    ir_[1].w() = 0.25;

    ir_[2].x() = pa;
    ir_[2].y() = pa;
    ir_[2].z() = 0.0;
    ir_[2].w() = 0.25;

    ir_[3].x() = na;
    ir_[3].y() = pa;
    ir_[3].z() = 0.0;
    ir_[3].w() = 0.25;
}


Foam::finiteElements::quadrilateral4Shell2::quadrilateral4Shell2()
{
    ir_.setSize(8);

    const scalar a = 1.0/sqrt(3.0);
    const scalar na = (1.0 - a)*0.5;
    const scalar pa = 1.0 - na;
    const scalar w = 0.125;

    ir_[0].x() = na;
    ir_[0].y() = na;
    ir_[0].z() = -a;
    ir_[0].w() = w;

    ir_[1].x() = pa;
    ir_[1].y() = na;
    ir_[1].z() = -a;
    ir_[1].w() = w;

    ir_[2].x() = pa;
    ir_[2].y() = pa;
    ir_[2].z() = -a;
    ir_[2].w() = w;

    ir_[3].x() = na;
    ir_[3].y() = pa;
    ir_[3].z() = -a;
    ir_[3].w() = w;

    ir_[4].x() = na;
    ir_[4].y() = na;
    ir_[4].z() = a;
    ir_[4].w() = w;

    ir_[5].x() = a;
    ir_[5].y() = na;
    ir_[5].z() = a;
    ir_[5].w() = w;

    ir_[6].x() = pa;
    ir_[6].y() = pa;
    ir_[6].z() = a;
    ir_[6].w() = w;

    ir_[7].x() = na;
    ir_[7].y() = pa;
    ir_[7].z() = a;
    ir_[7].w() = w;
}


Foam::finiteElements::quadrilateral4Shell3::quadrilateral4Shell3()
{
    ir_.setSize(12);

    const scalar a = 1.0/sqrt(3.0);
    const scalar na = (1.0 - a)*0.5;
    const scalar pa = 1.0 - na;
    const scalar w1 = 5.0/45.0;
    const scalar w2 = 1.0/9.0;

    ir_[0].x() = na;
    ir_[0].y() = na;
    ir_[0].z() = -a;
    ir_[0].w() = w1;

    ir_[1].x() = pa;
    ir_[1].y() = na;
    ir_[1].z() = -a;
    ir_[1].w() = w1;

    ir_[2].x() = pa;
    ir_[2].y() = pa;
    ir_[2].z() = -a;
    ir_[2].w() = w1;

    ir_[3].x() = na;
    ir_[3].y() = pa;
    ir_[3].z() = -a;
    ir_[3].w() = w1;

    ir_[4].x() = na;
    ir_[4].y() = na;
    ir_[4].z() = 0.0;
    ir_[4].w() = w2;

    ir_[5].x() = a;
    ir_[5].y() = na;
    ir_[5].z() = 0.0;
    ir_[5].w() = w2;

    ir_[6].x() = pa;
    ir_[6].y() = pa;
    ir_[6].z() = 0.0;
    ir_[6].w() = w2;

    ir_[7].x() = na;
    ir_[7].y() = pa;
    ir_[7].z() = 0.0;
    ir_[7].w() = w2;

    ir_[8].x() = na;
    ir_[8].y() = na;
    ir_[8].z() = a;
    ir_[8].w() = w1;

    ir_[9].x() = pa;
    ir_[9].y() = na;
    ir_[9].z() = a;
    ir_[9].w() = w1;

    ir_[10].x() = pa;
    ir_[10].y() = pa;
    ir_[10].z() = a;
    ir_[10].w() = w1;

    ir_[11].x() = na;
    ir_[11].y() = pa;
    ir_[11].z() = a;
    ir_[11].w() = w1;
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::finiteElements::quadrilateral4Shell::~quadrilateral4Shell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::finiteElements::quadrilateral4Shell::getNodes
(
    const List<vector>& verts
) const
{
    return tmp<pointField>(new pointField(verts));
}


Foam::scalarList Foam::finiteElements::quadrilateral4Shell::calcShape
(
    const vector& pt
) const
{
    scalarList shape(this->nNodes());
    shape[0] = pt.x()*pt.y();
    shape[1] = (1.0 - pt.x())*pt.y();
    shape[2] = (1.0 - pt.x())*(1.0 - pt.y());
    shape[3] = pt.x()*(1.0 - pt.y());

    return shape;
}


Foam::scalarRectangularMatrix Foam::finiteElements::quadrilateral4Shell::calcDShape
(
    const vector& pt
) const
{
    scalarRectangularMatrix dshape(this->nNodes(), this->nDims());
    dshape[0][0] =  pt.y();
    dshape[1][0] = -pt.y();
    dshape[2][0] = -(1.0 - pt.y());
    dshape[3][0] = (1.0 - pt.y());

    dshape[0][1] =  pt.x();
    dshape[1][1] =  (1.0 - pt.x());
    dshape[2][1] = -(1.0 - pt.x());
    dshape[3][1] = -pt.x();

    return dshape;
}


void Foam::finiteElements::quadrilateral4Shell::vtkData
(
    labelList& data,
    label& start,
    const labelList& labels
) const
{
    forAll(labels, i)
    {
        data[start++] = labels[i];
    }
}


// ************************************************************************* //
