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

#include "feMesh.H"
#include "Time.H"
#include "SubList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::feMesh::calcElements() const
{
    if (debug)
    {
        InfoInFunction << "Assembling shape functions" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (elementsPtr_)
    {
        FatalErrorInFunction
            << "elements already exist"
            << abort(FatalError);
    }

    elementsPtr_ = new List<element>(mesh().nCells());
    List<element>& elements = *elementsPtr_;
    forAll(elements, elemi)
    {
        elements[elemi].initialize(mesh(), elemi, order_, GeoType::CELL);
    }

    if (order_ > 1)
    {
        if (nodesPtr_ || edgeNodesPtr_ || faceNodesPtr_)
        {
            FatalErrorInFunction
                << "Nodes already exist"
                << abort(FatalError);
        }

        // References to poly mesh data
        const pointField& points = mesh().points();
        const labelListList& cellEdges = mesh().cellEdges();
        const cellList& cells = mesh().cells();

        // Dynamic list of nodes
        DynamicList<vector> nodes(points);

        //- Nodes associated with edges and faces
        edgeNodesPtr_ = new List<labelList>(mesh().nEdges());
        faceNodesPtr_ = new List<labelList>(mesh().nFaces());
        List<labelList>& edgeNodes = *edgeNodesPtr_;
        List<labelList>& faceNodes = *faceNodesPtr_;

        // Temporary storage
        DynamicList<label> newElem;
        DynamicList<point> elemNodes;

        forAll(elements, celli)
        {
            element& elem = elements[celli];
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
                const labelList& cEdges = cellEdges[celli];
                label edgei = -1;
                forAll(cEdges, ej)
                {
                    if (e == edges[cEdges[ej]])
                    {
                        edgei = cEdges[ej];
                        break;
                    }
                }
                if (edgei == -1)
                {
                    FatalErrorInFunction
                        << "Could not find edge " << e << " in mesh" << endl
                        << abort(FatalError);
                }

                // Insert edges nodes if not already created
                if (!edgeNodes[edgei].size())
                {
                    edgeNodes[edgei].setSize(nEdgeNodes);
                    nodes.append(SubList<point>(elemNodes, nEdgeNodes, ni));
                    for (label eni = 0; eni < nEdgeNodes; eni++)
                    {
                        newElem.append(nodes.size() + eni);
                        edgeNodes[edgei][eni] = nodes.size() + eni;
                    }
                }
                else
                {
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
                }

                // Increment starting index
                ni = newElem.size();
            }

            // Add face nodes
            const faceList faces(elem.faces());
            forAll(faces, fi)
            {
                const face& f = faces[fi];

                // Find the corresponding polyMesh face
                const cell& cFaces = cells[celli];
                label facei = -1;
                forAll(cFaces, fj)
                {
                    if (f == faces[cFaces[fj]])
                    {
                        facei = cFaces[fj];
                        break;
                    }
                }
                if (facei == -1)
                {
                    FatalErrorInFunction
                        << "Could not find face " << f << " in mesh" << endl
                        << abort(FatalError);
                }

                // Number of internal face nodes for the given face
                const label nFaceNodes = elem.nFaceNodes(fi);

                // Insert face nodes if not already exisiting
                if (!faceNodes[facei].size())
                {
                    faceNodes[facei].setSize(nFaceNodes);
                    nodes.append(SubList<point>(elemNodes, nFaceNodes, ni));
                    for (label fni = 0; fni < nFaceNodes; fni++)
                    {
                        newElem.append(nodes.size());
                        faceNodes[facei][fni] = nodes.size();
                    }
                }
                else
                {
                    // Add nodes already existing by matching
                    newElem.append
                    (
                        finiteElement::matchNodes
                        (
                            SubList<vector>(elemNodes, nFaceNodes, ni),
                            nodes,
                            faceNodes[facei]
                        )
                    );
                }

                // Incremenet the node index
                ni = newElem.size();
            }


            // Add internal nodes
            for (ni = newElem.size(); ni < elemNodes.size(); ni++)
            {
                newElem.append(nodes.size());
                nodes.append(elemNodes[ni]);
            }

            // Transfer the node indices to the element
            elem.transfer(newElem);
        }
        nodesPtr_ = new pointIOField
        (
            IOobject
            (
                "nodes",
                objectRegistry::time().timeName(),
                polyMesh::meshSubDir,
                mesh()
            ),
            pointField(move(nodes))
        );
    }
}


void Foam::feMesh::calcShapes() const
{
    if (debug)
    {
        InfoInFunction << "Assembling shape functions" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (shapesPtr_)
    {
        FatalErrorInFunction
            << "shape functions areas already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    shapesPtr_ = new List<List<scalarList>>(elements.size());
    List<List<scalarList>>& shapes = *shapesPtr_;
    forAll(shapes, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        shapes[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            shapes[elemi][i] = curElem.calcShape(ir[i]);
        }
    }
}


void Foam::feMesh::calcDShapes() const
{
    if (debug)
    {
        InfoInFunction << "Assembling dshape functions" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (dshapesPtr_)
    {
        FatalErrorInFunction
            << "shape function derivatives areas already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    dshapesPtr_ = new List<List<scalarRectangularMatrix>>(elements.size());
    List<List<scalarRectangularMatrix>>& dshapes = *dshapesPtr_;
    forAll(dshapes, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        dshapes[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            dshapes[elemi][i] = curElem.calcDShape(ir[i]);
        }
    }
}


void Foam::feMesh::calcBs() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling physical shape function gradients" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (BsPtr_)
    {
        FatalErrorInFunction
            << "shape function derivatives areas already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<scalarRectangularMatrix>>& ds = this->dshapes();
    const List<List<tensor>>& iJ = invJs();
    BsPtr_ = new List<List<scalarRectangularMatrix>>(elements.size());
    List<List<scalarRectangularMatrix>>& Bs = *BsPtr_;
    forAll(*this, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        Bs[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            scalarRectangularMatrix& curB = Bs[elemi][i];
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
}


void Foam::feMesh::calcJs() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling jacobians" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (JsPtr_)
    {
        FatalErrorInFunction
            << "Jacobains already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const pointField& nodes = this->nodes();
    JsPtr_ = new List<List<tensor>>(elements.size());
    List<List<tensor>>& Js = *JsPtr_;
    forAll(Js, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        Js[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            Js[elemi][i] = curElem.calcJ(ir[i], nodes);
        }
    }
}


void Foam::feMesh::calcInvJs() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling inverse jacobians" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (invJsPtr_)
    {
        FatalErrorInFunction
            << "Inverse jacobains already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<tensor>>& Js = this->Js();

    invJsPtr_ = new List<List<tensor>>(elements.size());
    List<List<tensor>>& invJs = *invJsPtr_;
    forAll(invJs, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        invJs[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            invJs[elemi][i] = curElem.calcInvJ(Js[elemi][i]);
        }
    }
}


void Foam::feMesh::calcWs() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling element weights" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (WsPtr_)
    {
        FatalErrorInFunction
            << "Element weights already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<tensor>>& Js = this->Js();

    WsPtr_ = new List<List<scalar>>(elements.size());
    List<List<scalar>>& Ws = *WsPtr_;
    forAll(Ws, elemi)
    {
        const element& curElem = elements[elemi];
        const integrationRule& ir = curElem.ir();
        Ws[elemi].setSize(ir.size());
        forAll(ir, i)
        {
            Ws[elemi][i] = curElem.calcW(Js[elemi][i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::element>& Foam::feMesh::elements() const
{
    if (!elementsPtr_)
    {
        calcElements();
    }

    return *elementsPtr_;
}


const Foam::List<Foam::labelList>& Foam::feMesh::edgeNodes() const
{
    if (!edgeNodesPtr_)
    {
        calcElements();
    }

    return *edgeNodesPtr_;
}


const Foam::List<Foam::labelList>& Foam::feMesh::faceNodes() const
{
    if (!faceNodesPtr_)
    {
        calcElements();
    }

    return *faceNodesPtr_;
}


const Foam::pointField& Foam::feMesh::nodes() const
{
    if (order_ < 2)
    {
        return mesh().points();
    }

    if (!nodesPtr_)
    {
        calcElements();
    }

    return *nodesPtr_;
}


const Foam::List<Foam::List<Foam::scalarList>>& Foam::feMesh::shapes() const
{
    if (!shapesPtr_)
    {
        calcShapes();
    }

    return *shapesPtr_;
}


const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::feMesh::dshapes() const
{
    if (!dshapesPtr_)
    {
        calcDShapes();
    }

    return *dshapesPtr_;
}


const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::feMesh::Bs() const
{
    if (!BsPtr_)
    {
        calcBs();
    }

    return *BsPtr_;
}


const Foam::List<Foam::List<Foam::tensor>>& Foam::feMesh::Js() const
{
    if (!JsPtr_)
    {
        calcJs();
    }

    return *JsPtr_;
}

const Foam::List<Foam::List<Foam::tensor>>& Foam::feMesh::invJs() const
{
    if (!invJsPtr_)
    {
        calcInvJs();
    }

    return *invJsPtr_;
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::feMesh::Ws() const
{
    if (!WsPtr_)
    {
        calcWs();
    }

    return *WsPtr_;
}


// ************************************************************************* //
