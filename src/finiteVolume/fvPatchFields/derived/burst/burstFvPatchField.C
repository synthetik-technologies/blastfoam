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

#include "burstFvPatchField.H"
#include "burstFvPatchBase.H"
#include "coupledFvPatchField.H"
#include "zeroGradientFvPatchField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstFvPatchField<Type>::burstFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    burstFvPatchFieldBase(p),
    burstPatchField_
    (
        new calculatedFvPatchField<Type>
        (
            p,
            iF
        )
    ),
    intactPatchField_
    (
        new calculatedFvPatchField<Type>
        (
            p,
            iF
        )
    )
{
    if (notNull(iF))
    {
        this->operator=(this->patchInternalField());
    }
}


template<class Type>
Foam::burstFvPatchField<Type>::burstFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    burstFvPatchFieldBase(p),
    burstPatchField_(),
    intactPatchField_()
{
    if (!isA<burstFvPatchBase>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "    patch type '" << p.type()
            << "' not derived from '" << burstFvPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalIOError);
    }

    // Create a new patch dictionary and replace the type with the intactType
    {
        dictionary burstDict(dict.parent(), dict.subDict("burstPatch"));
        if (!burstDict.found("patchType"))
        {
            burstDict.add("patchType", this->patch().type());
        }
        burstPatchField_.set
        (
            fvPatchField<Type>::New
            (
                p,
                iF,
                burstDict
            ).ptr()
        );
    }
    {
        dictionary intactDict(dict.parent(), dict.subDict("intactPatch"));
        if (!intactDict.found("patchType"))
        {
            intactDict.add("patchType", this->patch().type());
        }
        intactPatchField_.set
        (
            fvPatchField<Type>::New
            (
                p,
                iF,
                intactDict
            ).ptr()
        );
    }

    Field<Type>::operator=
    (
        intact()*intactPatchField_()
      + (1.0 - intact())*burstPatchField_()
    );
}


template<class Type>
Foam::burstFvPatchField<Type>::burstFvPatchField
(
    const burstFvPatchField<Type>& bpf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(bpf, p, iF, mapper),
    burstFvPatchFieldBase(p, bpf, mapper),
    burstPatchField_
    (
        fvPatchField<Type>::New
        (
            bpf.burstPatchField_(),
            p,
            iF,
            mapper
        ).ptr()
    ),
    intactPatchField_
    (
        fvPatchField<Type>::New
        (
            bpf.intactPatchField_(),
            p,
            iF,
            mapper
        ).ptr()
    )
{}


template<class Type>
Foam::burstFvPatchField<Type>::burstFvPatchField
(
    const burstFvPatchField<Type>& bpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(bpf, iF),
    burstFvPatchFieldBase(this->patch()),
    burstPatchField_(bpf.burstPatchField_->clone(iF).ptr()),
    intactPatchField_(bpf.intactPatchField_->clone(iF).ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::patchNeighbourField() const
{
    return
        intact()
       *(
            intactPatchField_->coupled()
          ? intactPatchField_().patchNeighbourField()
          : intactPatchField_()
        )
      + (1.0 - intact())
       *(
            burstPatchField_->coupled()
          ? burstPatchField_().patchNeighbourField()
          : burstPatchField_()
        );
}


template<class Type>
void Foam::burstFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    intactPatchField_->autoMap(m);
    burstPatchField_->autoMap(m);
}


template<class Type>
void Foam::burstFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const burstFvPatchField<Type>& bpf =
        refCast<const burstFvPatchField<Type>>(ptf);
    intactPatchField_->rmap(bpf.intactPatchField_(), addr);
    burstPatchField_->rmap(bpf.burstPatchField_(), addr);
}


template<class Type>
void Foam::burstFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Before we call evaluate on the cyclic patch, assign the current field
    // to the intact patch. This is done to make sure everything that has been
    // done to the actual patch is transfered
    intactPatchField_() = *this;
    burstPatchField_() = *this;

    intactPatchField_->updateCoeffs();
    burstPatchField_->updateCoeffs();

    burstFvPatchFieldBase::update();

    fvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::burstFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    burstPatchField_->initEvaluate(commsType);
    intactPatchField_->initEvaluate(commsType);

    fvPatchField<Type>::initEvaluate(commsType);
}


template<class Type>
void Foam::burstFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    intactPatchField_->evaluate(commsType);
    burstPatchField_->evaluate(commsType);

    Field<Type>::operator=
    (
        intact()*intactPatchField_()
      + (1.0 - intact())*burstPatchField_()
    );

    fvPatchField<Type>::evaluate(commsType);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_->valueInternalCoeffs(w)*intact()
      + burstPatchField_->valueInternalCoeffs(w)*intact();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_->valueBoundaryCoeffs(w)*intact()
      + burstPatchField_->valueBoundaryCoeffs(w)*(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        (
            intactPatchField_->coupled()
          ? intactPatchField_->gradientInternalCoeffs(deltaCoeffs)
          : intactPatchField_->gradientInternalCoeffs()
        )*intact()
      + (
            burstPatchField_->coupled()
          ? burstPatchField_->gradientInternalCoeffs(deltaCoeffs)
          : burstPatchField_->gradientInternalCoeffs()
        )*(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientInternalCoeffs() const
{
    return
        (
            intactPatchField_->coupled()
          ? intactPatchField_->gradientInternalCoeffs
            (
                this->patch().deltaCoeffs()
            )
          : intactPatchField_->gradientInternalCoeffs()
        )*intact()
      + (
            burstPatchField_->coupled()
          ? burstPatchField_->gradientInternalCoeffs
            (
                this->patch().deltaCoeffs()
            )
          : burstPatchField_->gradientInternalCoeffs()
        )*(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        (
            intactPatchField_->coupled()
          ? intactPatchField_->gradientBoundaryCoeffs(deltaCoeffs)
          : intactPatchField_->gradientBoundaryCoeffs()
        )*intact()
      + (
            burstPatchField_->coupled()
          ? burstPatchField_->gradientBoundaryCoeffs(deltaCoeffs)
          : burstPatchField_->gradientBoundaryCoeffs()
        )*(1.0 - intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        (
            intactPatchField_->coupled()
          ? intactPatchField_->gradientBoundaryCoeffs
            (
                this->patch().deltaCoeffs()
            )
          : intactPatchField_->gradientBoundaryCoeffs()
        )*intact()
      + (
            burstPatchField_->coupled()
          ? burstPatchField_->gradientBoundaryCoeffs
            (
                this->patch().deltaCoeffs()
            )
          : burstPatchField_->gradientBoundaryCoeffs()
        )*(1.0 - intact());
}


template<class Type>
void Foam::burstFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    if (isA<coupledFvPatchField<Type>>(burstPatchField_()))
    {
        scalarField bResult(result.size(), Zero);
        refCast<const coupledFvPatchField<Type>>
        (
            burstPatchField_()
        ).updateInterfaceMatrix
        (
            bResult,
            psiInternal,
            coeffs,
            cmpt,
            commsType
        );
        result -= (1.0 - intact())*bResult;
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
        result -= iResult*intact();
    }
}


template<class Type>
void Foam::burstFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    if (isA<coupledFvPatchField<Type>>(burstPatchField_()))
    {
        Field<Type> bResult(result.size(), Zero);
        refCast<const coupledFvPatchField<Type>>
        (
            burstPatchField_()
        ).updateInterfaceMatrix
        (
            bResult,
            psiInternal,
            coeffs,
            commsType
        );
        result -= (1.0 - intact())*bResult;
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
        result -= intact()*iResult;
    }
}


template<class Type>
void Foam::burstFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    {
        // Writing is a little weird since the burstPatchField has a different
        // type, but is in the same dictionary
        OStringStream oss;
        burstPatchField_->write(oss);
        dictionary dict(IStringStream(oss.str())());
        os.indent();
        os << "burstPatch" << dict;
    }
    {
        // Writing is a little weird since the intactPatchField has a different
        // type, but is in the same dictionary
        OStringStream oss;
        intactPatchField_->write(oss);
        dictionary dict(IStringStream(oss.str())());
        os.indent();
        os << "intactPatch" << dict;
    }
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    fvPatchField<Type>::operator=(ul);
    intactPatchField_() = ul;
    burstPatchField_() = ul;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=(ptf);
    intactPatchField_() = ptf;
    burstPatchField_() = ptf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator+=(ptf);
    intactPatchField_() += ptf;
    burstPatchField_() += ptf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator-=(ptf);
    intactPatchField_() -= ptf;
    burstPatchField_() -= ptf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<Type>::operator*=(ptf);
    intactPatchField_() *= ptf;
    burstPatchField_() *= ptf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<Type>::operator/=(ptf);
    intactPatchField_() /= ptf;
    burstPatchField_() /= ptf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
    intactPatchField_() += tf;
    burstPatchField_() += tf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
    intactPatchField_() -= tf;
    burstPatchField_() -= tf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
    intactPatchField_() *= tf;
    burstPatchField_() *= tf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
    intactPatchField_() /= tf;
    burstPatchField_() /= tf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_() = t;
    burstPatchField_() = t;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
    intactPatchField_() += t;
    burstPatchField_() += t;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
    intactPatchField_() -= t;
    burstPatchField_() -= t;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
    intactPatchField_() *= s;
    burstPatchField_() *= s;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
    intactPatchField_() /= s;
    burstPatchField_() /= s;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
    intactPatchField_() == ptf;
    burstPatchField_() == ptf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
    intactPatchField_() == tf;
    burstPatchField_() == tf;
}


template<class Type>
void Foam::burstFvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_() == t;
    burstPatchField_() == t;
}

// ************************************************************************* //
