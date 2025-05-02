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

#include "triangleShellFiniteElements.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace finiteElements
{
    defineTypeNameAndDebug(triangle3Shell1, 0);
    defineTypeNameAndDebug(triangle3Shell2, 0);
    defineTypeNameAndDebug(triangle3Shell3, 0);

    addToRunTimeSelectionTable(shellFiniteElement, triangle3Shell1, type);
    addToRunTimeSelectionTable(shellFiniteElement, triangle3Shell2, type);
    addToRunTimeSelectionTable(shellFiniteElement, triangle3Shell3, type);


    Foam::List<Foam::integrationRule>
        Foam::triangle3Shell::integrationRules;
}
}


// // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
//
// Foam::finiteElements::triangle3Shell::triangle3Shell()
// :
//     FiniteElement<ElementType::TRI>(1, 3)
// {
//     nodes_[0][0] = 0.0;
//     nodes_[0][1] = 0.0;
//
//     nodes_[1][0] = 1.0;
//     nodes_[1][1] = 0.0;
//
//     nodes_[2][0] = 0.0;
//     nodes_[2][1] = 1.0;
// }
//
//
// Foam::finiteElements::triangle3Shell1::triangle3Shell1()
// {
//     ir_.setSize(3);
//
//     const scalar a = 1.0/6.0;
//     const scalar b = 2.0/3.0;
//     const scalar w = 1.0/3.0;
//
//     ir_[0].x() = a;
//     ir_[0].y() = a;
//     ir_[0].z() = 0.0;
//     ir_[0].w() = w;
//
//     ir_[1].x() = b;
//     ir_[1].y() = a;
//     ir_[1].z() = 0.0;
//     ir_[1].w() = w;
//
//     ir_[2].x() = a;
//     ir_[2].y() = b;
//     ir_[2].z() = 0.0;
//     ir_[2].w() = w;
// }
//
//
// Foam::finiteElements::triangle3Shell2::triangle3Shell2()
// {
//     ir_.setSize(6);
//
//     const scalar a = 1.0/6.0;
//     const scalar b = 2.0/3.0;
//     const scalar c = 1.0/sqrt(3.0);
//     const scalar w = 1.0/6.0;
//
//     ir_[0].x() = a;
//     ir_[0].y() = a;
//     ir_[0].z() = -c;
//     ir_[0].w() = w;
//
//     ir_[1].x() = b;
//     ir_[1].y() = a;
//     ir_[1].z() = -c;
//     ir_[1].w() = w;
//
//     ir_[2].x() = a;
//     ir_[2].y() = b;
//     ir_[2].z() = -c;
//     ir_[2].w() = w;
//
//     ir_[3].x() = a;
//     ir_[3].y() = a;
//     ir_[3].z() = c;
//     ir_[3].w() = w;
//
//     ir_[4].x() = b;
//     ir_[4].y() = a;
//     ir_[4].z() = c;
//     ir_[4].w() = w;
//
//     ir_[5].x() = a;
//     ir_[5].y() = b;
//     ir_[5].z() = c;
//     ir_[5].w() = w;
// }
//
//
// Foam::finiteElements::triangle3Shell3::triangle3Shell3()
// {
//     ir_.setSize(6);
//
//     const scalar a = 1.0/6.0;
//     const scalar b = 2.0/3.0;
//     const scalar w1 = 5.0/9.0;
//     const scalar w2 = 8.0/9.0;
//
//     ir_[0].x() = a;
//     ir_[0].y() = a;
//     ir_[0].z() = -b;
//     ir_[0].w() = a*w1;
//
//     ir_[1].x() = b;
//     ir_[1].y() = a;
//     ir_[1].z() = -b;
//     ir_[1].w() = a*w1;
//
//     ir_[2].x() = a;
//     ir_[2].y() = b;
//     ir_[2].z() = -b;
//     ir_[2].w() = a*w1;
//
//     ir_[3].x() = a;
//     ir_[3].y() = a;
//     ir_[3].z() = 0.0;
//     ir_[3].w() = a*w2;
//
//     ir_[4].x() = b;
//     ir_[4].y() = a;
//     ir_[4].z() = 0.0;
//     ir_[4].w() = a*w2;
//
//     ir_[5].x() = a;
//     ir_[5].y() = b;
//     ir_[5].z() = 0.0;
//     ir_[5].w() = a*w2;
//
//     ir_[6].x() = a;
//     ir_[6].y() = a;
//     ir_[6].z() = b;
//     ir_[6].w() = a*w1;
//
//     ir_[7].x() = b;
//     ir_[7].y() = a;
//     ir_[7].z() = b;
//     ir_[7].w() = a*w1;
//
//     ir_[8].x() = a;
//     ir_[8].y() = b;
//     ir_[8].z() = b;
//     ir_[8].w() = a*w1;
// }
//
// // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
//
// Foam::finiteElements::triangle3Shell::~triangle3Shell()
// {}
//
// // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
//
// Foam::tmp<Foam::pointField> Foam::finiteElements::triangle3Shell::getNodes
// (
//     const vector& p0,
//     const vector& p1,
//     const vector& p2
// ) const
// {
//     tmp<pointField> tnodes(new pointField(this->nNodes()));
//     pointField& nodes = tnodes.ref();
//
//     nodes[0] = p0;
//     nodes[1] = p1;
//     nodes[2] = p2;
//
//     return tnodes;
// }
//
//
// Foam::tmp<Foam::pointField> Foam::finiteElements::triangle3Shell::getNodes
// (
//     const List<vector>& verts
// ) const
// {
//     return getNodes(verts[0], verts[1], verts[2]);
// }
//
//
// Foam::scalarList Foam::finiteElements::triangle3Shell::calcShape
// (
//     const vector& pt
// ) const
// {
//     if (order_ == 0)
//     {
//         return scalarList(this->nNodes(), 1.0);
//     }
//
//     Field<scalar> u(this->nNodes());
//     u[0] = 1.0 - pt.x() - pt.y();
//     u[1] = pt.x();
//     u[2] = pt.y();
//     return u;
// }
//
//
// Foam::scalarRectangularMatrix Foam::finiteElements::triangle3Shell::calcDShape
// (
//     const vector& pt
// ) const
// {
//     scalarRectangularMatrix du(this->nNodes(), this->nDims());
//     du[0][0] = -1.0;
//     du[1][0] = 1.0;
//     du[2][0] = 0.0;
//
//     du[0][1] = -1.0;
//     du[1][1] = 0.0;
//     du[2][2] = 0.0;
//     return du;
// }
//
//
// void Foam::finiteElements::triangle3Shell::calcDShape
// (
//     const vector& pt,
//     scalarList& shape,
//     scalarRectangularMatrix& dshape
// ) const
// {
//
//     shape[0] = 1.0 - pt.x() - pt.y();
//     shape[1] = pt.x();
//     shape[2] = pt.y();
//
//     dshape[0][0] = -1.0;
//     dshape[1][0] = 1.0;
//     dshape[2][0] = 0.0;
//
//     dshape[0][1] = -1.0;
//     dshape[1][1] = 0.0;
//     dshape[2][2] = 0.0;
// }
//
//
// void Foam::finiteElements::triangle3Shell::vtkData
// (
//     labelList& data,
//     label& start,
//     const labelList& labels
// ) const
// {
//     forAll(labels, i)
//     {
//         data[start++] = labels[i];
//     }
// }
//
//
// //- Return integration rule for a given order
// const Foam::integrationRule& Foam::triangle3Shell::ir
// (
//     const label order
// ) const
// {
//     if (integrationRules.size() <= order)
//     {
//         integrationRules.setSize(order+1);
//     }
//
//     integrationRule& ir = integrationRules[order];
//     if (!ir.size())
//     {
//         // Triangle integration
//         const integrationRules& ir2 =
//             integrationRule::triRule(2);
//
//         // Integration through thickness
//         const integrationRules& irs =
//             integrationRule::segRule(order);
//
//         const label n = ir2.size()*irs.size();
//         ir.setSize(n);
//         label I = 0;
//         forAll(ir2, i)
//         {
//             forAll(irs, j)
//             {
//                 ir[I++].set3
//                 (
//                     ir2[i].x(),
//                     ir2[i].y(),
//                     2.0*irs[j].x() - 1.0;
//                     ir2[i].w()*irs[j].w()
//                 );
//             }
//         }
//     }
//     return ir;

// ************************************************************************* //
