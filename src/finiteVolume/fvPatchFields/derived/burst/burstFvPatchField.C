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
    fvPatchField<Type>(p, iF),
    burstPatchField_
    (
        new zeroGradientFvPatchField<Type>
        (
            p,
            iF
        )
    ),
    intactPatchField_
    (
        new zeroGradientFvPatchField<Type>
        (
            p,
            iF
        )
    ),
    pName_("p"),
    impulseName_("impulse"),
    burstPatch_
    (
        const_cast<burstFvPatch&>
        (
            dynamicCast<const burstFvPatch>(p)
        )
    )
{
    Field<Type>::operator=(Zero);
    burstPatchField_() = Zero;
    intactPatchField_() = Zero;
}


template<class Type>
Foam::burstFvPatchField<Type>::burstFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    burstPatchField_(),
    intactPatchField_(),
    pName_(dict.lookupOrDefault<word>("pName", "p")),
    impulseName_(dict.lookupOrDefault<word>("impulseName", "impulse")),
    burstPatch_
    (
        const_cast<burstFvPatch&>
        (
            refCast<const burstFvPatch>(p)
        )
    )
{
    if (!isA<burstFvPatch>(p))
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
        if (!burstPatch_.uptoDate(this->db().time().timeIndex()))
        {
            burstPatch_.intact() =
                scalarField("intact", dict, p.size());
        }
    }
    if (dict.found("pRef"))
    {
        burstPatch_.pRef() =
            scalarField("pRef", dict, p.size());
    }

    // Create a new patch dictionary and replace the type with the intactType
    {
        dictionary burstDict(dict.parent(), dict.subDict("burstPatch"));
        if (!burstDict.found("patchType"))
        {
            burstDict.add("patchType", typeName);
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
    }

    Field<Type>::operator=
    (
        burstPatch_.intact()*intactPatchField_()
      + (1.0 - burstPatch_.intact())*burstPatchField_()
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
    fvPatchField<Type>(bpf, p, iF, mapper),
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
    ),
    pName_(bpf.pName_),
    impulseName_(bpf.impulseName_),
    burstPatch_
    (
        const_cast<burstFvPatch&>
        (
            refCast<const burstFvPatch>(p)
        )
    )
{}


template<class Type>
Foam::burstFvPatchField<Type>::burstFvPatchField
(
    const burstFvPatchField<Type>& bpf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(bpf, iF),
    burstPatchField_(bpf.burstPatchField_->clone(iF).ptr()),
    intactPatchField_(bpf.intactPatchField_->clone(iF).ptr()),
    pName_(bpf.pName_),
    impulseName_(bpf.impulseName_),
    burstPatch_
    (
        const_cast<burstFvPatch&>
        (
            refCast<const burstFvPatch>(this->patch())
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::patchNeighbourField() const
{
    const scalarField& intact = burstPatch_.intact();
    return
        intact*intactPatchField_().patchNeighbourField()
      + (1.0 - intact)*burstPatchField_().patchNeighbourField();
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
    burstPatch_.autoMap(m);
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
    burstPatch_.rmap(ptf.patch(), addr);
}


template<class Type>
void Foam::burstFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    intactPatchField_->updateCoeffs();
    burstPatchField_->updateCoeffs();

    // Execute the change to the openFraction only once per time-step
    if (!burstPatch_.uptoDate(this->db().time().timeIndex()))
    {
        scalarField p(this->patch().patch().size(), 0.0);
        scalarField impulse(this->patch().patch().size(), 0.0);

        const polyMesh& mesh = this->patch().boundaryMesh().mesh();
        bool shouldUpdate = false;

        if (burstPatch_.usePressure())
        {
            if (mesh.template foundObject<volScalarField>(pName_))
            {
                const burstFvPatchField<scalar>& pp =
                    refCast<const burstFvPatchField<scalar>>
                    (
                        this->patch().
                        template lookupPatchField<volScalarField, scalar>
                        (
                            pName_
                        )
                    );
                p = mag(pp.patchInternalField());
                shouldUpdate = true;
            }
            else
            {
                WarningInFunction
                    << "Could not find " << pName_
                    << ", neglecting pressure. " << endl;
            }
        }
        if (burstPatch_.useImpulse())
        {
            if (mesh.template foundObject<volScalarField>(impulseName_))
            {
                const burstFvPatchField<scalar>& pImpulse =
                    refCast<const burstFvPatchField<scalar>>
                    (
                        this->patch().
                        template lookupPatchField<volScalarField, scalar>
                        (
                            impulseName_
                        )
                    );
                impulse = mag(pImpulse.patchInternalField());
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
            burstPatch_.update(p, impulse);
        }
    }
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

    // Before we call evaluate on the cyclic patch, assign the current field
    // to the intact patch. This is done to make sure everything that has been
    // done to the actual patch is transfered
    intactPatchField_() = *this;
    burstPatchField_() = *this;

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
        burstPatch_.intact()*intactPatchField_()
      + (1.0 - burstPatch_.intact())*burstPatchField_()
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
        intactPatchField_->valueInternalCoeffs(w)
       *burstPatch_.intact()
      + burstPatchField_->valueInternalCoeffs(w)
       *burstPatch_.intact();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return
        intactPatchField_->valueBoundaryCoeffs(w)
       *burstPatch_.intact()
      + burstPatchField_->valueBoundaryCoeffs(w)
       *(1.0 - burstPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        intactPatchField_->gradientInternalCoeffs(deltaCoeffs)
       *burstPatch_.intact()
      + burstPatchField_->gradientInternalCoeffs(deltaCoeffs)
       *(1.0 - burstPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientInternalCoeffs() const
{
    return
        intactPatchField_->gradientInternalCoeffs()
       *burstPatch_.intact()
      + burstPatchField_->gradientInternalCoeffs()
       *(1.0 - burstPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return
        intactPatchField_->gradientBoundaryCoeffs(deltaCoeffs)
       *burstPatch_.intact()
      + burstPatchField_->gradientBoundaryCoeffs(deltaCoeffs)
       *(1.0 - burstPatch_.intact());
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        intactPatchField_->gradientBoundaryCoeffs()
       *burstPatch_.intact()
      + burstPatchField_->gradientBoundaryCoeffs()
       *(1.0 - burstPatch_.intact());
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
    const scalarField& intact = burstPatch_.intact();
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
        result -= (1.0 - intact)*bResult;
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
void Foam::burstFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    const scalarField& intact = burstPatch_.intact();
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
        result -= (1.0 - intact)*bResult;
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
void Foam::burstFvPatchField<Type>::write(Ostream& os) const
{
    os  << "burstPatch" << nl
        << token::BEGIN_BLOCK << nl
        << incrIndent << burstPatchField_() << decrIndent
        << token::END_BLOCK;
    os  << "intactPatch" << nl
        << token::BEGIN_BLOCK << nl
        << incrIndent << intactPatchField_() << decrIndent
        << token::END_BLOCK;

    fvPatchField<Type>::write(os);
    writeEntryIfDifferent<word>(os, "pName", "p", pName_);
    writeEntryIfDifferent<word>(os, "impulseName", "impulse", impulseName_);

    writeEntry(os, "intact", burstPatch_.intact());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

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
