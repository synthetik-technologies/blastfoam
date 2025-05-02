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

#include "feMesh1.H"
#include "Time.H"
#include "SubList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::feMesh1::calcElements() const
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

    if (shell_)
    {
        elementsPtr_ = new List<element>(mesh().nFaces());
        List<element>& elements = *elementsPtr_;
        forAll(elements, elemi)
        {
            elements[elemi].initializeShell
            (
                mesh(),
                elemi,
                1, // Element order
                intOrder_, // Integration order
                shellIntOrder_,
                GeoType::FACE
            );
        }
    }
    else
    {
        elementsPtr_ = new List<element>(mesh().nCells());
        List<element>& elements = *elementsPtr_;
        forAll(elements, elemi)
        {
            elements[elemi].initialize
            (
                mesh(),
                elemi,
                1, // Element order
                intOrder_, // Integration order
                GeoType::CELL
            );
        }
    }
}


void Foam::feMesh1::calcIpLabels() const
{
    if (debug)
    {
        InfoInFunction << "Assembling intrgration point labels" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ipLabelsPtr_)
    {
        FatalErrorInFunction
            << "Integration point labels already exist"
            << abort(FatalError);
    }

    const List<element>& elements = this->elements();;

    ipLabelsPtr_ = new List<labelList>(elements.size());
    List<labelList>& ipLabels = *ipLabelsPtr_;
    nIp_ = 0;
    forAll(elements, elemi)
    {
        const label nip = elements[elemi].ir().size();
        ipLabels[elemi].setSize(nip);
        for (label i = 0; i < nip; i++)
        {
            ipLabels[elemi][i] = nIp_++;
        }
    }

    // forAll(boundary_, patchi)
    // {
    //     nIp_ += boundary_[patchi].updateIpLabels(nIp_);
    // }
}


void Foam::feMesh1::calcShapes() const
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


void Foam::feMesh1::calcDShapes() const
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


void Foam::feMesh1::calcBs() const
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
    forAll(elements, elemi)
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

            for (label si = 0; si < dshape.m(); si++)
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


void Foam::feMesh1::calcJs() const
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


void Foam::feMesh1::calcInvJs() const
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


void Foam::feMesh1::calcWs() const
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


void Foam::feMesh1::calcW() const
{
    if (debug)
    {
        InfoInFunction
            << "Assembling element weight matrix" << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (WPtr_)
    {
        FatalErrorInFunction
            << "Element weights already exist"
            << abort(FatalError);
    }

    const elementList& elements = this->elements();
    const List<List<scalarList>>& shapes = this->shapes();
    const List<List<scalar>>& Ws = this->Ws();

    WPtr_ = new pointScalarField
    (
        IOobject
        (
            "W",
            mesh().time().name(),
            mesh()
        ),
        pointMesh_,
        0.0
    );
    pointScalarField& W = *WPtr_;

    forAll(elements, ei)
    {
        const element& elem = elements[ei];

        UIndirectList<scalar> w_loc(W, elem);

        const integrationRule& ir = elem.ir();
        forAll(ir, rulei)
        {
            const integrationPoint& ip = ir[rulei];

            const scalarList& shape = shapes[ei][rulei];
            const scalar w = ip.w()*Ws[ei][rulei];

            forAll(shape, si)
            {
                forAll(shape, sj)
                {
                    w_loc[si] += w*shape[si]*shape[sj];
                }
            }
        }
    }
    syncTools::syncPointList(mesh(), W, plusEqOp<scalar>(), 0.0);

    if (gMin(W.primitiveField()) < small)
    {
        FatalErrorInFunction
            << "Weight matrix is singular" << endl
            << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::List<Foam::element>& Foam::feMesh1::elements() const
{
    if (!elementsPtr_)
    {
        calcElements();
    }

    return *elementsPtr_;
}


Foam::label Foam::feMesh1::nIp() const
{
    if (!ipLabelsPtr_)
    {
        calcIpLabels();
    }

    return nIp_;
}


const Foam::List<Foam::labelList>& Foam::feMesh1::ipLabels() const
{
    if (!ipLabelsPtr_)
    {
        calcIpLabels();
    }

    return *ipLabelsPtr_;
}


const Foam::pointField& Foam::feMesh1::nodes() const
{
    return mesh().points();
}


const Foam::List<Foam::List<Foam::scalarList>>& Foam::feMesh1::shapes() const
{
    if (!shapesPtr_)
    {
        calcShapes();
    }

    return *shapesPtr_;
}


const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::feMesh1::dshapes() const
{
    if (!dshapesPtr_)
    {
        calcDShapes();
    }

    return *dshapesPtr_;
}


const Foam::List<Foam::List<Foam::scalarRectangularMatrix>>&
Foam::feMesh1::Bs() const
{
    if (!BsPtr_)
    {
        calcBs();
    }

    return *BsPtr_;
}


const Foam::List<Foam::List<Foam::tensor>>& Foam::feMesh1::Js() const
{
    if (!JsPtr_)
    {
        calcJs();
    }

    return *JsPtr_;
}

const Foam::List<Foam::List<Foam::tensor>>& Foam::feMesh1::invJs() const
{
    if (!invJsPtr_)
    {
        calcInvJs();
    }

    return *invJsPtr_;
}

const Foam::List<Foam::List<Foam::scalar>>& Foam::feMesh1::Ws() const
{
    if (!WsPtr_)
    {
        calcWs();
    }

    return *WsPtr_;
}


const Foam::pointScalarField& Foam::feMesh1::W() const
{
    if (!WPtr_)
    {
        calcW();
    }

    return *WPtr_;
}

// ************************************************************************* //
