/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "calculatedProcessorFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::calculatedProcessorFvPatchField<Type>::calculatedProcessorFvPatchField
(
    const lduInterface& interface,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    //lduInterfaceField(interface),
    coupledFvPatchField<Type>(p, iF),
    procInterface_(refCast<const lduPrimitiveProcessorInterface>(interface)),
    sendBuf_(interface.faceCells().size()),
    receiveBuf_(interface.faceCells().size()),
    scalarSendBuf_(interface.faceCells().size()),
    scalarReceiveBuf_(interface.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


template<class Type>
Foam::calculatedProcessorFvPatchField<Type>::calculatedProcessorFvPatchField
(
    const calculatedProcessorFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    //lduInterfaceField(ptf.procInterface_),
    coupledFvPatchField<Type>(ptf, iF),
    procInterface_(ptf.procInterface_),
    sendBuf_(procInterface_.faceCells().size()),
    receiveBuf_(procInterface_.faceCells().size()),
    scalarSendBuf_(procInterface_.faceCells().size()),
    scalarReceiveBuf_(procInterface_.faceCells().size()),
    outstandingSendRequest_(-1),
    outstandingRecvRequest_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::calculatedProcessorFvPatchField<Type>::ready() const
{
    if
    (
        this->outstandingSendRequest_ >= 0
     && this->outstandingSendRequest_ < Pstream::nRequests()
    )
    {
        bool finished =
            UPstream::finishedRequest(this->outstandingSendRequest_);
        if (!finished)
        {
            return false;
        }
    }
    this->outstandingSendRequest_ = -1;

    if
    (
        this->outstandingRecvRequest_ >= 0
     && this->outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        bool finished =
            UPstream::finishedRequest(this->outstandingRecvRequest_);
        if (!finished)
        {
            return false;
        }
    }
    this->outstandingRecvRequest_ = -1;

    return true;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::calculatedProcessorFvPatchField<Type>::patchNeighbourField() const
{
    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "On patch of size " << procInterface_.faceCells().size()
            << " between proc " << procInterface_.myProcNo()
            << " and " << procInterface_.neighbProcNo()
            << " outstanding request."
            << abort(FatalError);
    }
    return *this;
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        //this->patchInternalField(sendBuf_);
        // Bypass patchInternalField since uses fvPatch addressing
        {
            const Field<Type>& iF = this->internalField();
            const labelList& fc = procInterface_.faceCells();
            sendBuf_.setSize(fc.size());
            forAll(fc, i)
            {
                sendBuf_[i] = iF[fc[i]];
            }
        }

        // Receive straight into *this
        this->setSize(sendBuf_.size());
        outstandingRecvRequest_ = UPstream::nRequests();
        UIPstream::read
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<char*>(this->data()),
            this->byteSize(),
            procInterface_.tag(),
            procInterface_.comm()
        );

        outstandingSendRequest_ = UPstream::nRequests();
        UOPstream::write
        (
            Pstream::commsTypes::nonBlocking,
            procInterface_.neighbProcNo(),
            reinterpret_cast<const char*>(sendBuf_.cdata()),
            this->byteSize(),
            procInterface_.tag(),
            procInterface_.comm()
        );
    }
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (Pstream::parRun())
    {
        if
        (
            outstandingRecvRequest_ >= 0
         && outstandingRecvRequest_ < Pstream::nRequests()
        )
        {
            UPstream::waitRequest(outstandingRecvRequest_);
        }
        outstandingSendRequest_ = -1;
        outstandingRecvRequest_ = -1;
    }
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::initInterfaceMatrixUpdate
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    // Bypass patchInternalField since uses fvPatch addressing
    const labelList& fc = procInterface_.faceCells();
    scalarSendBuf_.setSize(fc.size());
    forAll(fc, i)
    {
        scalarSendBuf_[i] = psiInternal[fc[i]];
    }

    if (debug && !this->ready())
    {
        FatalErrorInFunction
            << "On patch " //<< interface_.name()
            << " outstanding request."
            << abort(FatalError);
    }


    scalarReceiveBuf_.setSize(scalarSendBuf_.size());
    outstandingRecvRequest_ = UPstream::nRequests();
    UIPstream::read
    (
        Pstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        reinterpret_cast<char*>(scalarReceiveBuf_.data()),
        scalarReceiveBuf_.byteSize(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    outstandingSendRequest_ = UPstream::nRequests();
    UOPstream::write
    (
        Pstream::commsTypes::nonBlocking,
        procInterface_.neighbProcNo(),
        reinterpret_cast<const char*>(scalarSendBuf_.cdata()),
        scalarSendBuf_.byteSize(),
        procInterface_.tag(),
        procInterface_.comm()
    );

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = false;
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::addToInternalField
(
    scalarField& result,
    const bool add,
    const scalarField& coeffs,
    const scalarField& vals
) const
{
    const labelUList& faceCells = this->procInterface_.faceCells();

    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*vals[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*vals[elemI];
        }
    }
}


template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
//     const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        UPstream::waitRequest(outstandingRecvRequest_);
    }
    // Recv finished so assume sending finished as well.
    outstandingSendRequest_ = -1;
    outstandingRecvRequest_ = -1;

    // Consume straight from scalarReceiveBuf_. Note use of our own
    // helper to avoid using fvPatch addressing
//     addToInternalField(result, !add, coeffs, scalarReceiveBuf_);
    addToInternalField(result, true, coeffs, scalarReceiveBuf_);

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = true;
}

template<class Type>
void Foam::calculatedProcessorFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    if (this->updatedMatrix())
    {
        return;
    }

    if
    (
        outstandingRecvRequest_ >= 0
     && outstandingRecvRequest_ < Pstream::nRequests()
    )
    {
        UPstream::waitRequest(outstandingRecvRequest_);
    }
    // Recv finished so assume sending finished as well.
    outstandingSendRequest_ = -1;
    outstandingRecvRequest_ = -1;

    // Consume straight from scalarReceiveBuf_. Note use of our own
    // helper to avoid using fvPatch addressing
//     addToInternalField(result, !add, coeffs, scalarReceiveBuf_);
//     addToInternalField(result, true, coeffs, scalarReceiveBuf_);

    const_cast<lduInterfaceField&>
    (
        static_cast<const lduInterfaceField&>(*this)
    ).updatedMatrix() = true;
}


//template<class Type>
//void Foam::calculatedProcessorFvPatchField<Type>::write(Ostream& os) const
//{
//    //zeroGradientFvPatchField<Type>::write(os);
//    fvPatchField<Type>::write(os);
//    this->writeEntry("value", os);
//}


// ************************************************************************* //
