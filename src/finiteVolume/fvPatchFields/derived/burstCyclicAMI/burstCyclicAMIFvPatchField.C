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

#include "burstCyclicAMIFvPatchField.H"
#include "zeroGradientFvPatchField.H"
#include "volFields.H"
#include "burstCyclicAMIFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstCyclicAMIFvPatchField<Type>::burstCyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMIFvPatchField<Type>(p, iF),
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
Foam::burstCyclicAMIFvPatchField<Type>::burstCyclicAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicAMIFvPatchField<Type>(p, iF, dict),
    burstFvPatchFieldBase(p),
    intactPatchField_()
{
    if (!isA<burstCyclicAMIFvPatch>(p))
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
}


template<class Type>
Foam::burstCyclicAMIFvPatchField<Type>::burstCyclicAMIFvPatchField
(
    const burstCyclicAMIFvPatchField<Type>& bpf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicAMIFvPatchField<Type>(bpf, p, iF, mapper),
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
Foam::burstCyclicAMIFvPatchField<Type>::burstCyclicAMIFvPatchField
(
    const burstCyclicAMIFvPatchField<Type>& bpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicAMIFvPatchField<Type>(bpf, iF),
    burstFvPatchFieldBase(this->patch()),
    intactPatchField_(bpf.intactPatchField_->clone(iF))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicAMIFvPatchField<Type>::patchNeighbourField() const
{
    tmp<Field<Type>> pnF(cyclicAMIFvPatchField<Type>::patchNeighbourField());
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
    scalarField in(intact());
    return in*inF + (1.0 - in)*pnF;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    cyclicAMIFvPatchField<Type>::autoMap(m);
    intactPatchField_->autoMap(m);
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    cyclicAMIFvPatchField<Type>::rmap(ptf, addr);

    const burstCyclicAMIFvPatchField<Type>& bpf =
        refCast<const burstCyclicAMIFvPatchField<Type>>(ptf);
    intactPatchField_->rmap(bpf.intactPatchField_(), addr);
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Set the intact field to zero gradient in the event that the
    // patchField has not been updated
    intactPatchField_.ref() = this->patchInternalField();

    intactPatchField_->updateCoeffs();
    cyclicAMIFvPatchField<Type>::updateCoeffs();
    burstFvPatchFieldBase::update();
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    cyclicAMIFvPatchField<Type>::initEvaluate(commsType);
    intactPatchField_->initEvaluate(commsType);
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    intactPatchField_->evaluate(commsType);
    if (block_)
    {
        Field<Type>::operator=(intactPatchField_());
    }
    else
    {
        cyclicAMIFvPatchField<Type>::evaluate(commsType);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicAMIFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_().valueInternalCoeffs(w)*intact()
      + cyclicAMIFvPatchField<Type>::valueInternalCoeffs(w)
       *(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicAMIFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_().valueBoundaryCoeffs(w)*intact()
      + cyclicAMIFvPatchField<Type>::valueBoundaryCoeffs(w)
       *(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicAMIFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        (
            intactPatchField_().coupled()
          ? intactPatchField_().gradientInternalCoeffs(deltaCoeffs)
          : intactPatchField_().gradientInternalCoeffs()
        )*intact()
      + cyclicAMIFvPatchField<Type>::gradientInternalCoeffs(deltaCoeffs)
       *(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstCyclicAMIFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        (
            intactPatchField_().coupled()
          ? intactPatchField_().gradientBoundaryCoeffs(deltaCoeffs)
          : intactPatchField_().gradientBoundaryCoeffs()
        )*intact()
      + cyclicAMIFvPatchField<Type>::gradientBoundaryCoeffs(deltaCoeffs)
       *(1.0 - intact());
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
        cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
void Foam::burstCyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
        cyclicAMIFvPatchField<Type>::updateInterfaceMatrix
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
void Foam::burstCyclicAMIFvPatchField<Type>::write(Ostream& os) const
{
    cyclicAMIFvPatchField<Type>::write(os);
    {
        // Writing is a little weird since the intactPatchField has a different
        // type, but is in the same dictionary
        OStringStream oss;
        intactPatchField_->write(oss);
        dictionary dict(IStringStream(oss.str())());
        os.indent();
        os << "intactPatch" << dict;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    Field<Type>::operator=(ul);
    intactPatchField_.ref() = ul;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
    intactPatchField_.ref() = ptf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator+=(ptf);
    intactPatchField_.ref() += ptf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator-=(ptf);
    intactPatchField_.ref() -= ptf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    Field<Type>::operator*=(ptf);
    intactPatchField_.ref() *= ptf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    Field<Type>::operator/=(ptf);
    intactPatchField_.ref() /= ptf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
    intactPatchField_.ref() += tf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
    intactPatchField_.ref() -= tf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
    intactPatchField_.ref() *= tf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
    intactPatchField_.ref() /= tf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_.ref() = t;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
    intactPatchField_.ref() += t;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
    intactPatchField_.ref() -= t;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
    intactPatchField_.ref() *= s;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
    intactPatchField_.ref() /= s;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
    intactPatchField_.ref() == ptf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
    intactPatchField_.ref() == tf;
}


template<class Type>
void Foam::burstCyclicAMIFvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_.ref() == t;
}

// ************************************************************************* //
