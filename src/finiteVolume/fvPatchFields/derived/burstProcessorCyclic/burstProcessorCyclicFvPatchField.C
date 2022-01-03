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

#include "burstProcessorCyclicFvPatchField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstProcessorCyclicFvPatchField<Type>::burstProcessorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(p, iF),
    intactPatchField_
    (
        fvPatchField<Type>::New
        (
            "zeroGradient",
            typeName,
            p,
            iF
        ).ptr()
    ),
    pName_("p"),
    impulseName_("impulse"),
    burstProcessorCyclicPatch_
    (
        const_cast<burstProcessorCyclicFvPatch&>
        (
            dynamicCast<const burstProcessorCyclicFvPatch>(p)
        )
    )
{
    this->operator=(Zero);
}


template<class Type>
Foam::burstProcessorCyclicFvPatchField<Type>::burstProcessorCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchField<Type>(p, iF, dict),
    intactPatchField_(),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    burstProcessorCyclicPatch_
    (
        const_cast<burstProcessorCyclicFvPatch&>
        (
            refCast<const burstProcessorCyclicFvPatch>(p)
        )
    )
{
    if (!isA<burstProcessorCyclicFvPatch>(p))
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

    if (dict.found("intact"))
    {
        if (!burstProcessorCyclicPatch_.uptoDate(this->db().time().timeIndex()))
        {
            burstProcessorCyclicPatch_.intact() =
                scalarField("intact", dict, p.size());
        }
    }

    // Create a new patch dictionary and replace the type with the intactType
    dictionary intactDict(dict.parent(), dict);
    intactDict.set("type", dict.lookup<word>("intactType"));
    if (!intactDict.found("patchType"))
    {
        intactDict.add("patchType", typeName);
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

    this->evaluate(Pstream::commsTypes::blocking);
}


template<class Type>
Foam::burstProcessorCyclicFvPatchField<Type>::burstProcessorCyclicFvPatchField
(
    const burstProcessorCyclicFvPatchField<Type>& bpf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvPatchField<Type>(bpf, p, iF, mapper),
    intactPatchField_(bpf.intactPatchField_->clone(iF).ptr()),
    pName_(bpf.pName_),
    impulseName_(bpf.impulseName_),
    burstProcessorCyclicPatch_
    (
        const_cast<burstProcessorCyclicFvPatch&>
        (
            refCast<const burstProcessorCyclicFvPatch>(p)
        )
    )
{}


template<class Type>
Foam::burstProcessorCyclicFvPatchField<Type>::burstProcessorCyclicFvPatchField
(
    const burstProcessorCyclicFvPatchField<Type>& bpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(bpf, iF),
    intactPatchField_(bpf.intactPatchField_->clone(iF).ptr()),
    pName_(bpf.pName_),
    impulseName_(bpf.impulseName_),
    burstProcessorCyclicPatch_
    (
        const_cast<burstProcessorCyclicFvPatch&>
        (
            refCast<const burstProcessorCyclicFvPatch>(this->patch())
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::patchNeighbourField() const
{
    const scalarField& intact = burstProcessorCyclicPatch_.intact();
    return
        intact*intactPatchField_()
      + (1.0 - intact)*cyclicFvPatchField<Type>::patchNeighbourField();
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    cyclicFvPatchField<Type>::autoMap(m);
    intactPatchField_->autoMap(m);
    burstProcessorCyclicPatch_.autoMap(m);
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    cyclicFvPatchField<Type>::rmap(ptf, addr);

    const burstProcessorCyclicFvPatchField<Type>& bpf =
        refCast<const burstProcessorCyclicFvPatchField<Type>>(ptf);
    intactPatchField_->rmap(bpf.intactPatchField_(), addr);
    burstProcessorCyclicPatch_.rmap(ptf.patch(), addr);
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    intactPatchField_->updateCoeffs();
    cyclicFvPatchField<Type>::updateCoeffs();

    // Execute the change to the openFraction only once per time-step
    if (!burstProcessorCyclicPatch_.uptoDate(this->db().time().timeIndex()))
    {
        scalarField p(this->patch().patch().size(), 0.0);
        scalarField impulse(this->patch().patch().size(), 0.0);

        const polyMesh& mesh = this->patch().boundaryMesh().mesh();
        bool shouldUpdate = false;

        if (burstProcessorCyclicPatch_.usePressure())
        {
            if (mesh.template foundObject<volScalarField>(pName_))
            {
                const burstProcessorCyclicFvPatchField<scalar>& pp =
                    refCast<const burstProcessorCyclicFvPatchField<scalar>>
                    (
                        this->patch().
                        template lookupPatchField<volScalarField, scalar>
                        (
                            pName_
                        )
                    );
                p =
                    mag
                    (
                        pp.cyclicFvPatchField<scalar>::patchNeighbourField()
                      - pp.patchInternalField()
                    );
                shouldUpdate = true;
            }
            else
            {
                WarningInFunction
                    << "Could not find " << pName_
                    << ", neglecting pressure. " << endl;
            }
        }
        if (burstProcessorCyclicPatch_.useImpulse())
        {
            if (mesh.template foundObject<volScalarField>(impulseName_))
            {
                const burstProcessorCyclicFvPatchField<scalar>& pImpulse =
                    refCast<const burstProcessorCyclicFvPatchField<scalar>>
                    (
                        this->patch().
                        template lookupPatchField<volScalarField, scalar>
                        (
                            impulseName_
                        )
                    );
                impulse =
                    mag
                    (
                        pImpulse.cyclicFvPatchField<scalar>::patchNeighbourField()
                      - pImpulse.patchInternalField()
                    );
                shouldUpdate = true;
            }
            else
            {
                WarningInFunction
                    << "Could not find " << impulseName_
                    << ", neglecting impulse. " << endl;
            }
        }

        if (shouldUpdate)
        {
            burstProcessorCyclicPatch_.update(p, impulse);
        }
    }
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    cyclicFvPatchField<Type>::initEvaluate(commsType);
    intactPatchField_->initEvaluate(commsType);
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Before we call evaluate on the cyclic patch, assign the current field
    // to the intact patch. This is done to make sure everything that has been
    // done to the actual patch is transfered
    intactPatchField_() = *this;

    intactPatchField_->evaluate(commsType);
    cyclicFvPatchField<Type>::evaluate(commsType);

    Field<Type>::operator=
    (
        burstProcessorCyclicPatch_.intact()*intactPatchField_()
      + (1.0 - burstProcessorCyclicPatch_.intact())*(*this)
    );

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_().valueInternalCoeffs(w)
       *burstProcessorCyclicPatch_.intact()
      + cyclicFvPatchField<Type>::valueInternalCoeffs(w)
       *burstProcessorCyclicPatch_.intact();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_().valueBoundaryCoeffs(w)
       *burstProcessorCyclicPatch_.intact()
      + cyclicFvPatchField<Type>::valueBoundaryCoeffs(w)
       *(1.0 - burstProcessorCyclicPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        intactPatchField_().gradientInternalCoeffs(deltaCoeffs)
       *burstProcessorCyclicPatch_.intact()
      + cyclicFvPatchField<Type>::gradientInternalCoeffs(deltaCoeffs)
       *(1.0 - burstProcessorCyclicPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::gradientInternalCoeffs() const
{
    return
        intactPatchField_().gradientInternalCoeffs()
       *burstProcessorCyclicPatch_.intact()
      + cyclicFvPatchField<Type>::gradientInternalCoeffs()
       *(1.0 - burstProcessorCyclicPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        intactPatchField_().gradientBoundaryCoeffs(deltaCoeffs)
       *burstProcessorCyclicPatch_.intact()
      + cyclicFvPatchField<Type>::gradientBoundaryCoeffs(deltaCoeffs)
       *(1.0 - burstProcessorCyclicPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstProcessorCyclicFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        intactPatchField_().gradientBoundaryCoeffs()
       *burstProcessorCyclicPatch_.intact()
      + cyclicFvPatchField<Type>::gradientBoundaryCoeffs()
       *(1.0 - burstProcessorCyclicPatch_.intact());
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    const scalarField& intact = burstProcessorCyclicPatch_.intact();
    {
        scalarField cResult(result.size(), Zero);
        cyclicFvPatchField<Type>::updateInterfaceMatrix
        (
            cResult,
            psiInternal,
            coeffs,
            cmpt,
            commsType
        );
        result -= (1.0 - intact)*cResult;
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
        result -= iResult*intact;
    }
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    const scalarField& intact = burstProcessorCyclicPatch_.intact();
    {
        Field<Type> cResult(result.size(), Zero);
        cyclicFvPatchField<Type>::updateInterfaceMatrix
        (
            cResult,
            psiInternal,
            coeffs,
            commsType
        );
        result -= (1.0 - intact)*cResult;
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
        result -= intact*iResult;
    }
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::write(Ostream& os) const
{
    {
        // Writing is a little weird since the intactPatchField has a different
        // type, but is in the same dictionary
        OStringStream oss;
        intactPatchField_->write(oss);
        dictionary dict(IStringStream(oss.str())());

        dict.changeKeyword("type", "intactType", false);
        dict.remove("patchType");
        forAllConstIter(IDLList<entry>, dict, iter)
        {
            iter().write(os);
        }
    }

    fvPatchField<Type>::write(os);
    writeEntryIfDifferent<word>(os, "pName", "p", pName_);
    writeEntryIfDifferent<word>(os, "impulseName", "impulse", impulseName_);

    writeEntry(os, "intact", burstProcessorCyclicPatch_.intact());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator=
(
    const UList<Type>& ul
)
{
    fvPatchField<Type>::operator=(ul);
    intactPatchField_() = ul;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=(ptf);
    intactPatchField_() = ptf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator+=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator+=(ptf);
    intactPatchField_() += ptf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator-=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator-=(ptf);
    intactPatchField_() -= ptf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator*=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<Type>::operator*=(ptf);
    intactPatchField_() *= ptf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator/=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<Type>::operator/=(ptf);
    intactPatchField_() /= ptf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator+=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator+=(tf);
    intactPatchField_() += tf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator-=
(
    const Field<Type>& tf
)
{
    Field<Type>::operator-=(tf);
    intactPatchField_() -= tf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator*=
(
    const scalarField& tf
)
{
    Field<Type>::operator*=(tf);
    intactPatchField_() *= tf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator/=
(
    const scalarField& tf
)
{
    Field<Type>::operator/=(tf);
    intactPatchField_() /= tf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator=
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_() = t;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator+=
(
    const Type& t
)
{
    Field<Type>::operator+=(t);
    intactPatchField_() += t;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator-=
(
    const Type& t
)
{
    Field<Type>::operator-=(t);
    intactPatchField_() -= t;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator*=
(
    const scalar s
)
{
    Field<Type>::operator*=(s);
    intactPatchField_() *= s;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator/=
(
    const scalar s
)
{
    Field<Type>::operator/=(s);
    intactPatchField_() /= s;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator==
(
    const fvPatchField<Type>& ptf
)
{
    Field<Type>::operator=(ptf);
    intactPatchField_() == ptf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator==
(
    const Field<Type>& tf
)
{
    Field<Type>::operator=(tf);
    intactPatchField_() == tf;
}


template<class Type>
void Foam::burstProcessorCyclicFvPatchField<Type>::operator==
(
    const Type& t
)
{
    Field<Type>::operator=(t);
    intactPatchField_() == t;
}

// ************************************************************************* //
