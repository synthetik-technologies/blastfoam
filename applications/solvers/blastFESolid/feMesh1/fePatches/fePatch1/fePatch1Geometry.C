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

#include "fePatch1.H"
#include "feMesh1.H"
#include "feBoundaryMesh1.H"
#include "Map.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fePatch1::calcElements() const
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

    const feMesh1& mesh = boundaryMesh_.mesh();
    const polyMesh& pmesh = mesh.mesh();

    elementsPtr_ = new List<element>(polyPatch_.size());
    List<element>& elements = *elementsPtr_;
    forAll(elements, elemi)
    {
        elements[elemi].initialize
        (
            pmesh,
            polyPatch_.start() + elemi,
            1, // Element order
            mesh.intOrder(), // Integration order
            GeoType::FACE
        );
    }
}


void Foam::fePatch1::calcLocalElements() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling local elements for patch "
            << polyPatch_.name() << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (localElementsPtr_)
    {
        FatalErrorInFunction
            << "localElements already exist"
            << abort(FatalError);
    }

    localElementsPtr_ = new List<labelList>(elements());
    List<labelList>& localElements = *localElementsPtr_;

    const Map<label>& meshPointMap = polyPatch_.meshPointMap();

    //- Local nodes
    forAll(localElements, ei)
    {
        // Create the local element for interpolation
        // on the patch
        labelList& localElem = localElements[ei];
        forAll(localElem, i)
        {
            localElem[i] = meshPointMap[localElem[i]];
        }
    }
}


void Foam::fePatch1::calcLocalNodes() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling local nodes for patch "
            << polyPatch_.name() << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (localNodesPtr_)
    {
        FatalErrorInFunction
            << "localNodes already exist"
            << abort(FatalError);
    }

    const Map<label>& meshPointMap = polyPatch_.meshPointMap();
    const pointField& nodes = boundaryMesh_.mesh().nodes();

    localNodesPtr_ = new pointField(meshPointMap.size());
    pointField& localNodes = *localNodesPtr_;

    //- Local nodes
    forAllConstIter(Map<label>, meshPointMap, iter)
    {
        localNodes[iter()] = nodes[iter.key()];
    }
}


void Foam::fePatch1::calcShapes() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcShapes() : "
               "calculating shape functions in fePatch1"
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
        Pout<< "fePatch1::calcShapes() : "
               "finished calculating shape functions in fePatch1"
            << endl;
    }
}


void Foam::fePatch1::calcDShapes() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcDShapes() : "
               "calculating derivatives of shape functions in fePatch1"
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
        Pout<< "fePatch1::calcDShapes() : "
               "finished calculating dertivatives of shape functions in fePatch1"
            << endl;
    }
}


void Foam::fePatch1::calcJs() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcJs() : "
               "calculating Jacobians in fePatch1"
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
        Pout<< "fePatch1::calcJs(): "
               "finished calculating Jacobians in fePatch1"
            << endl;
    }
}


void Foam::fePatch1::calcInvJs() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcInvJs(): "
               "calculating inverse Jacobians in fePatch1"
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
        Pout<< "fePatch1::calcInvJs(): "
               "finished calculating inverse Jacobians in fePatch1"
            << endl;
    }
}


void Foam::fePatch1::calcWs() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcWs(): "
               "calculating weights in fePatch1"
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
        Pout<< "fePatch1::calcWs(): "
               "finished calculating weights in fePatch1"
            << endl;
    }
}


void Foam::fePatch1::calcBs() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcBs() : "
               "calculating gradients of shape functions in fePatch1"
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
        Pout<< "fePatch1::calcBs() : "
               "finished calculating gradients of shape functions in fePatch1"
            << endl;
    }
}


void Foam::fePatch1::calcNodeNormals() const
{
    if (debug)
    {
        Pout<< "fePatch1::calcNodeNormals() : "
               "calculating nodeNormals in fePatch1"
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

    const elementList& elements = this->elements();
    const List<List<tensor>>& J = Js();

    nodeNormalsPtr_ = new List<List<vector>>(J.size());
    List<List<vector>>& ns = *nodeNormalsPtr_;

    forAll(J, elemi)
    {
        const element& curElem = elements[elemi];
        ns[elemi].setSize(J[elemi].size());

        forAll(ns[elemi], nodei)
        {
            ns[elemi][nodei] = curElem.calcOrtho(J[elemi][nodei]);
        }
    }

    if (debug)
    {
        Pout<< "fePatch1::calcNodeNormals() : "
               "finished calculating nodeNormals in fePatch1"
            << endl;
    }
}

// ************************************************************************* //
