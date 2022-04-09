/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "burstCyclicACMIFvPatchField.H"
#include "zeroGradientFvPatchField.H"
#include "volFields.H"
#include "burstCyclicACMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstCyclicACMIFvPatchField<Type>::burstCyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMIFvPatchField<Type>(p, iF),
    burstFvPatchFieldBase(p),
    intactPatchField_
    (
        new zeroGradientFvPatchField<Type>
        (
            p,
            iF
        )
    )
{
    this->operator=(Zero);
}


template<class Type>
Foam::burstCyclicACMIFvPatchField<Type>::burstCyclicACMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicACMIFvPatchField<Type>(p, iF, dict),
    burstFvPatchFieldBase(p),
    intactPatchField_()
{
    if (!isA<burstCyclicACMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not constraint type '" << typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    // Create a new patch dictionary and replace the type with the intactType
    dictionary intactDict(dict.parent(), dict.subDict("intactPatch"));
    if (!intactDict.found("patchType"))
    {
        intactDict.add("patchType", typeName);
    }
    intactPatchField_ =
        fvPatchField<Type>::New
        (
            p,
            iF,
            intactDict
        );

    cyclicACMIFvPatchField<Type>::evaluate(Pstream::commsTypes::blocking);
}


template<class Type>
Foam::burstCyclicACMIFvPatchField<Type>::burstCyclicACMIFvPatchField
(
    const burstCyclicACMIFvPatchField<Type>& bpf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicACMIFvPatchField<Type>(bpf, p, iF, mapper),
    burstFvPatchFieldBase(p),
    intactPatchField_
    (
        fvPatchField<Type>::New
        (
            bpf.intactPatchField_(),
            p,
            iF,
            mapper
        )
    )
{}


template<class Type>
Foam::burstCyclicACMIFvPatchField<Type>::burstCyclicACMIFvPatchField
(
    const burstCyclicACMIFvPatchField<Type>& bpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicACMIFvPatchField<Type>(bpf, iF),
    burstFvPatchFieldBase(this->patch()),
    intactPatchField_(bpf.intactPatchField_->clone(iF))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicACMIFvPatchField<Type>::patchNeighbourField() const
{
    tmp<Field<Type>> pnF(cyclicACMIFvPatchField<Type>::patchNeighbourField());
    tmp<Field<Type>> inF
    (
        intactPatchField_->coupled()
       ? intactPatchField_().patchNeighbourField()
       : intactPatchField_()
    );
    if (this->unblock_)
    {
        return pnF;
    }
    else if (this->block_)
    {
        return inF;
    }
    return intact()*inF + (1.0 - intact())*pnF;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    cyclicACMIFvPatchField<Type>::autoMap(m);
    intactPatchField_->autoMap(m);
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    cyclicACMIFvPatchField<Type>::rmap(ptf, addr);

    const burstCyclicACMIFvPatchField<Type>& bpf =
        refCast<const burstCyclicACMIFvPatchField<Type>>(ptf);
    intactPatchField_->rmap(bpf.intactPatchField_(), addr);
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Set the intact field to zero gradient in the event that the
    // patchField has not been updated
    intactPatchField_.ref() = this->patchInternalField();

    intactPatchField_->updateCoeffs();
    cyclicACMIFvPatchField<Type>::updateCoeffs();

    burstFvPatchFieldBase::update();
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    cyclicACMIFvPatchField<Type>::initEvaluate(commsType);
    intactPatchField_->initEvaluate(commsType);
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    intactPatchField_->evaluate(commsType);
    cyclicACMIFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicACMIFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_().valueInternalCoeffs(w)*intact()
      + cyclicACMIFvPatchField<Type>::valueInternalCoeffs(w)
       *intact();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicACMIFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_().valueBoundaryCoeffs(w)*intact()
      + cyclicACMIFvPatchField<Type>::valueBoundaryCoeffs(w)
       *(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicACMIFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        intactPatchField_().gradientInternalCoeffs()*intact()
      + cyclicACMIFvPatchField<Type>::gradientInternalCoeffs(deltaCoeffs)
       *(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicACMIFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        intactPatchField_().gradientBoundaryCoeffs()*intact()
      + cyclicACMIFvPatchField<Type>::gradientBoundaryCoeffs(deltaCoeffs)
       *(1.0 - intact());
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const polyPatch& p = this->patch().patch();
    const scalarField intact(this->intact());
    {
        scalarField cResult(result.size(), Zero);
        cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
        (
            cResult,
            psiInternal,
            coeffs,
            cmpt,
            commsType
        );
        forAll(p.faceCells(), fi)
        {
            const label celli = p.faceCells()[fi];
            result[celli] += (1.0 - intact[fi])*cResult[celli];
        }
    }

    if (isA<coupledFvPatchField<Type>>(intactPatchField_()))
    {
        scalarField iResult(result.size(), Zero);
        refCast<const coupledFvPatchField<Type>>
        (
            intactPatchField_()
        ).updateInterfaceMatrix
        (
            iResult,
            psiInternal,
            coeffs,
            cmpt,
            commsType
        );
        forAll(p.faceCells(), fi)
        {
            const label celli = p.faceCells()[fi];
            result[celli] += intact[fi]*iResult[celli];
        }
    }
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    const polyPatch& p = this->patch().patch();
    const scalarField intact(this->intact());
    {
        Field<Type> cResult(result.size(), Zero);
        cyclicACMIFvPatchField<Type>::updateInterfaceMatrix
        (
            cResult,
            psiInternal,
            coeffs,
            commsType
        );
        forAll(p.faceCells(), fi)
        {
            const label celli = p.faceCells()[fi];
            result[celli] += (1.0 - intact[fi])*cResult[celli];
        }
    }
    if (isA<coupledFvPatchField<Type>>(intactPatchField_()))
    {
        Field<Type> iResult(result.size(), Zero);
        refCast<const coupledFvPatchField<Type>>
        (
            intactPatchField_()
        ).updateInterfaceMatrix
        (
            iResult,
            psiInternal,
            coeffs,
            commsType
        );
        forAll(p.faceCells(), fi)
        {
            const label celli = p.faceCells()[fi];
            result[celli] += intact[fi]*iResult[celli];
        }
    }
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::write(Ostream& os) const
{
    cyclicACMIFvPatchField<Type>::write(os);

    writeKeyword(os, "intactPatch")
        << nl << indent << token::BEGIN_BLOCK << nl << incrIndent;
    intactPatchField_->write(os);
    os << decrIndent << indent << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
    intactPatchField_.ref() = ul;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
    intactPatchField_.ref() = ptf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator+=(ptf);
    intactPatchField_.ref() += ptf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator-=(ptf);
    intactPatchField_.ref() -= ptf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    Field<Type>::operator*=(ptf);
    intactPatchField_.ref() *= ptf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    Field<Type>::operator/=(ptf);
    intactPatchField_.ref() /= ptf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
    intactPatchField_.ref() += tf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
    intactPatchField_.ref() -= tf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
    intactPatchField_.ref() *= tf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
    intactPatchField_.ref() /= tf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_.ref() = t;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
    intactPatchField_.ref() += t;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
    intactPatchField_.ref() -= t;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
    intactPatchField_.ref() *= s;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
    intactPatchField_.ref() /= s;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
    intactPatchField_.ref() == ptf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
    intactPatchField_.ref() == tf;
}


template<class Type>
void Foam::burstCyclicACMIFvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_.ref() == t;
}

// ************************************************************************* //
