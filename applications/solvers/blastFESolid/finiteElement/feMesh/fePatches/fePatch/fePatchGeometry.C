/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "fePatch.H"
#include "feMesh.H"
#include "feBoundaryMesh.H"
#include "Map.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fePatch::calcElements() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling elements for patch "
            << polyPatch_.name() << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (elementsPtr_)
    {
        FatalErrorInFunction
            << "elements already exist"
            << abort(FatalError);
    }

    const feMesh& mesh = boundaryMesh_.mesh();
    const polyMesh& pmesh = mesh.mesh();
    const label order = mesh.order();

    elementsPtr_ = new List<element>(polyPatch_.size());
    List<element>& elements = *elementsPtr_;
    forAll(elements, elemi)
    {
        elements[elemi].initialize
        (
            pmesh,
            elemi,
            order,
            GeoType::FACE
        );
    }

    if (order > 1)
    {
        DynamicList<label> newElem;
        DynamicList<point> elemNodes;
        const labelListList& faceEdges = pmesh.faceEdges();
        const labelListList& edgeNodes = mesh.edgeNodes();
        const labelListList& faceNodes = mesh.faceNodes();
        const pointField& nodes = mesh.nodes();
        const pointField& points = pmesh.points();

        forAll(polyPatch_, fi)
        {
            const label facei = fi + polyPatch_.start();
            element& elem = elements[fi];
            newElem = elem;

            // Construct nodes for the given cell from the vertices
            elemNodes = elem.getNodes(pointField(points, elem));

            // Current node index
            label ni = newElem.size();

            // Number of edge nodes (internal)
            const label nEdgeNodes = elem.nEdgeNodes();

            // Get the edges of the element
            const edgeList edges(elem.edges());

            // Add all edge nodes
            forAll(edges, ei)
            {
                const edge& e = edges[ei];

                // Finde the corresponding mesh edge
                const labelList& fEdges = faceEdges[facei];                label edgei = -1;
                forAll(fEdges, ej)
                {
                    if (e == edges[fEdges[ej]])
                    {
                        edgei = fEdges[ej];
                        break;
                    }
                }
                if (edgei == -1)
                {
                    FatalErrorInFunction
                        << "Could not find edge " << e << " in mesh" << endl
                        << abort(FatalError);
                }

                // Add nodes already existing by matching
                newElem.append
                (
                    finiteElement::matchNodes
                    (
                        SubList<vector>(elemNodes, nEdgeNodes, ni),
                        nodes,
                        edgeNodes[edgei]
                    )
                );

                // Increment starting index
                ni = newElem.size();
            }

            // Add face nodes, elem is still the inital face
            newElem.append
            (
                finiteElement::matchNodes
                (
                    UIndirectList<point>(nodes, elem),
                    nodes,
                    faceNodes[facei]
                )
            );

            // Transfer indices
            elem.transfer(newElem);
        }
    }
}


void Foam::fePatch::calcMeshData() const
{
    if (debug)
    {
        Pout<< "fePatch::calcMeshData() : "
               "calculating mesh data in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate meshPoints
    // if they have already been calculated.
    if (meshNodesPtr_)
    {
        FatalErrorInFunction
            << "meshNodesPtr_ already allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();

    // Create a map for marking points.  Estimated size is 4 times the
    // number of faces in the patch
    labelHashSet markedNodes(4*elements.size());

    //- Unsorted version:
    DynamicList<label> meshNodes(2*elements.size());
    forAll(*this, elemi)
    {
        const labelList& curNodes = elements[elemi];

        forAll(curNodes, nodei)
        {
            if (markedNodes.insert(curNodes[nodei]))
            {
                meshNodes.append(curNodes[nodei]);
            }
        }
    }

    // Transfer to straight list (reuses storage)
    meshNodesPtr_ = new labelList(meshNodes, true);



    if (debug)
    {
        Pout<< "fePatch::calcMeshData() : "
               "finished calculating mesh data in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcShapes() const
{
    if (debug)
    {
        Pout<< "fePatch::calcShapes() : "
               "calculating shape functions in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (shapesPtr_)
    {
        FatalErrorInFunction
            << "shapesPts_ allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();

    shapesPtr_ = new List<List<scalarList>>(elements.size());
    List<List<scalarList>>& s = *shapesPtr_;

    forAll(elements, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        s[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            s[elemi][i] = curElem.calcShape(ir[i]);
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcShapes() : "
               "finished calculating shape functions in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcDShapes() const
{
    if (debug)
    {
        Pout<< "fePatch::calcDShapes() : "
               "calculating derivatives of shape functions in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (dshapesPtr_)
    {
        FatalErrorInFunction
            << "dshapesPtr_ allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();

    dshapesPtr_ = new List<List<scalarRectangularMatrix>>(elements.size());
    List<List<scalarRectangularMatrix>>& ds = *dshapesPtr_;

    forAll(elements, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        ds[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            ds[elemi][i] = curElem.calcDShape(ir[i]);
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcDShapes() : "
               "finished calculating dertivatives of shape functions in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcJs() const
{
    if (debug)
    {
        Pout<< "fePatch::calcJs() : "
               "calculating Jacobians in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (JsPtr_)
    {
        FatalErrorInFunction
            << "JPtr_ allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const pointField& nodes = boundaryMesh_.mesh().nodes();

    JsPtr_ = new List<List<tensor>>(elements.size());
    List<List<tensor>>& J = *JsPtr_;

    forAll(elements, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        J[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            J[elemi][i] = curElem.calcJ(ir[i], nodes);
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcJs(): "
               "finished calculating Jacobians in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcInvJs() const
{
    if (debug)
    {
        Pout<< "fePatch::calcInvJs(): "
               "calculating inverse Jacobians in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (invJsPtr_)
    {
        FatalErrorInFunction
            << "invJPtr_ allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<tensor>>& Js = this->Js();

    invJsPtr_ = new List<List<tensor>>(this->size());
    List<List<tensor>>& iJ = *invJsPtr_;

    forAll(elements, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        iJ[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            iJ[elemi][i] = curElem.calcInvJ(Js[elemi][i]);
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcInvJs(): "
               "finished calculating inverse Jacobians in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcWs() const
{
    if (debug)
    {
        Pout<< "fePatch::calcWs(): "
               "calculating weights in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (WsPtr_)
    {
        FatalErrorInFunction
            << "invJPtr_ allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<tensor>>& Js = this->Js();

    WsPtr_ = new List<List<scalar>>(elements.size());
    List<List<scalar>>& W = *WsPtr_;

    forAll(elements, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        W[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            W[elemi][i] = curElem.calcW(Js[elemi][i]);
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcWs(): "
               "finished calculating weights in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcBs() const
{
    if (debug)
    {
        Pout<< "fePatch::calcBs() : "
               "calculating gradients of shape functions in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate faceCentres
    // if they have already been calculated.
    if (BsPtr_)
    {
        FatalErrorInFunction
            << "BsPtr_ allocated"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<scalarRectangularMatrix>>& ds = dshapes();
    const List<List<tensor>>& iJ = invJs();

    BsPtr_ = new List<List<scalarRectangularMatrix>>(elements.size());
    List<List<scalarRectangularMatrix>>& B = *BsPtr_;

    forAll(elements, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        B[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            scalarRectangularMatrix& curB = B[elemi][i];
            const scalarRectangularMatrix& dshape = ds[elemi][i];
            const tensor& invJ = iJ[elemi][i];
            curB.setSize(dshape.m(), 3);
            curB = Zero;

            forAll(dshape, si)
            {
                for (label dimi = 0; dimi < dshape.n(); dimi++)
                {
                    for (label dimj = 0; dimj < 3; dimj++)
                    {
                        curB(si, dimj) += dshape(si, dimi)*invJ(dimi, dimj);
                    }
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcBs() : "
               "finished calculating gradients of shape functions in fePatch"
            << endl;
    }
}


void Foam::fePatch::calcNodeNormals() const
{
    if (debug)
    {
        Pout<< "fePatch::calcNodeNormals() : "
               "calculating nodeNormals in fePatch"
            << endl;
    }

    // It is considered an error to attempt to recalculate pointNormals
    // if they have already been calculated.
    if (nodeNormalsPtr_)
    {
        FatalErrorInFunction
            << "nodeNormalsPtr_already allocated"
            << abort(FatalError);
    }

    const List<List<tensor>>& J = Js();

    nodeNormalsPtr_ = new List<List<vector>>(J.size());
    List<List<vector>>& n = *nodeNormalsPtr_;

    forAll(J, elemi)
    {
        n[elemi].setSize(J[elemi].size());

        forAll(n[elemi], nodei)
        {
            const tensor& curJ = J[elemi][nodei];
            vector& curNormal = n[elemi][nodei];
            curNormal =
                vector
                (
                    curJ[1]*curJ[5] - curJ[2]*curJ[4],
                    curJ[2]*curJ[3] - curJ[0]*curJ[5],
                    curJ[0]*curJ[4] - curJ[1]*curJ[3]
                );
                curNormal /= mag(curNormal) + vSmall;
        }
    }

    if (debug)
    {
        Pout<< "fePatch::calcNodeNormals() : "
               "finished calculating nodeNormals in fePatch"
            << endl;
    }
}

// ************************************************************************* //
