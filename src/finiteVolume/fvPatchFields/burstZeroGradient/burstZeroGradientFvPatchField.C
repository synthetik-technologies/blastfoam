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

#include "burstZeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "cyclicFvPatch.H"
#include "cyclicAMIFvPatch.H"
#include "volFields.H"
#include "transformField.H"
#include "symmTensorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstZeroGradientFvPatchField<Type>::burstZeroGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    pName_("p"),
    cyclicPatchName_(),
    cyclicPatchLabel_(-1),
    pBurst_(great),
    burstImpulse_(-1.0),
    useImpulse_(false),
    impulse_(p.size(), 0.0),
    partialBurst_(false),
    intact_(p.size(), 1.0),
    curTimeIndex_(-1)
{}


template<class Type>
Foam::burstZeroGradientFvPatchField<Type>::burstZeroGradientFvPatchField
(
    const burstZeroGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    pBurst_(ptf.pBurst_),
    burstImpulse_(ptf.burstImpulse_),
    useImpulse_(ptf.useImpulse_),
    impulse_(mapper(ptf.impulse_)),
    partialBurst_(ptf.partialBurst_),
    intact_(mapper(ptf.intact_)),
    curTimeIndex_(ptf.curTimeIndex_)
{}


template<class Type>
Foam::burstZeroGradientFvPatchField<Type>::burstZeroGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    cyclicPatchName_(dict.lookup("cyclicPatch")),
    cyclicPatchLabel_(p.patch().boundaryMesh().findPatchID(cyclicPatchName_)),
    pBurst_(dict.lookup<scalar>("burstP")),
    burstImpulse_(dict.lookupOrDefault<scalar>("burstImpulse", -1)),
    useImpulse_(burstImpulse_ > 0),
    impulse_(p.size(), 0.0),
    partialBurst_(dict.lookupOrDefault("partialBurst", false)),
    intact_
    (
        dict.found("intact")
      ? scalarField("intact", dict, p.size())
      : scalarField(p.size(), 1.0)
    ),
    curTimeIndex_(-1)
{
    if (useImpulse_ && dict.found("impulse"))
    {
        impulse_ = scalarField("impulse", dict, p.size());
    }
    fvPatchField<Type>::operator=(this->patchInternalField());
}


template<class Type>
Foam::burstZeroGradientFvPatchField<Type>::burstZeroGradientFvPatchField
(
    const burstZeroGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    pName_(ptf.pName_),
    cyclicPatchName_(ptf.cyclicPatchName_),
    cyclicPatchLabel_(ptf.cyclicPatchLabel_),
    pBurst_(ptf.pBurst_),
    burstImpulse_(ptf.burstImpulse_),
    useImpulse_(ptf.useImpulse_),
    impulse_(ptf.impulse_),
    partialBurst_(ptf.partialBurst_),
    intact_(ptf.intact_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstZeroGradientFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    m(intact_, intact_);
    if (useImpulse_)
    {
        m(impulse_, impulse_);
    }

}


template<class Type>
void Foam::burstZeroGradientFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
    const burstZeroGradientFvPatchField<Type>& bbptf =
        refCast<const burstZeroGradientFvPatchField<Type>>(ptf);
    intact_.rmap(bbptf.intact_, addr);
    if (useImpulse_)
    {
        impulse_.rmap(bbptf.impulse_, addr);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::burstZeroGradientFvPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));

    if (isA<cyclicFvPatch>(this->patch()))
    {
        const cyclicFvPatch& cPatch =
            refCast<const cyclicFvPatch>(this->patch());
        const labelUList& nbrFaceCells =
            cPatch.cyclicPatch().nbrPatch().faceCells();
        Field<Type>& pnf = tpnf.ref();
        forAll(pnf, facei)
        {
            pnf[facei] = cPatch.transform().transform(iField[nbrFaceCells[facei]]);
        }
    }
    else if (isA<cyclicAMIFvPatch>(this->patch()))
    {
        const cyclicAMIFvPatch& cPatch =
            refCast<const cyclicAMIFvPatch>(this->patch());
        const labelUList& nbrFaceCells =
            cPatch.cyclicAMIPatch().nbrPatch().faceCells();

        Field<Type> pnf(iField, nbrFaceCells);
        if (cPatch.applyLowWeightCorrection())
        {
            tpnf =
                cPatch.interpolate(pnf, this->patchInternalField()());
        }
        else
        {
            tpnf = cPatch.interpolate(pnf);
        }
        cPatch.transform().transform(tpnf.ref(), tpnf());
    }

    return intact_*(*this) + (1.0 - intact_)*tpnf;
}


template<class Type>
void Foam::burstZeroGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    fixedValueFvPatchField<Type>::updateCoeffs();

    // Execute the change to the openFraction only once per time-step
    if
    (
        curTimeIndex_ != this->db().time().timeIndex()
     && this->db().template foundObject<volScalarField>(pName_)
    )
    {
        const volScalarField& p = this->db().template lookupObject<volScalarField>
        (
            pName_
        );
        scalarField pP
        (
            p.boundaryField()[this->patch().index()]
        );
        if (useImpulse_)
        {
            impulse_ += pP*this->db().time().deltaTValue();
        }

        if (partialBurst_)
        {
            intact_ = min(intact_, neg(pP - pBurst_));
            if (useImpulse_)
            {
                intact_ = min(intact_, neg(impulse_ - burstImpulse_));
            }
        }
        else
        {
            scalar maxP
            (
                max
                (
                    max(p.boundaryField()[this->patch().index()]),
                    max(p.boundaryField()[cyclicPatchLabel_])
                )
             );
            intact_ = min(intact_, scalar(neg(maxP - pBurst_)));
            if (useImpulse_)
            {
                intact_ =
                    min(intact_, scalar(neg(max(impulse_) - burstImpulse_)));
            }
        }

        typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
        FieldType& fld =
            this->db().template lookupObjectRef<FieldType>
            (
                this->internalField().name()
            );
        burstZeroGradientFvPatchField<Type>& nbr =
            dynamicCast<burstZeroGradientFvPatchField<Type>>
            (
                fld.boundaryFieldRef()[cyclicPatchLabel_]
            );
        nbr.updateCoeffs();
        intact_ = min(intact_, nbr.intact_);
        nbr.intact_ = intact_;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    tmp<vectorField> nHat = this->patch().nf();
    const Field<Type> iF(this->patchInternalField());
    Field<Type> intVals
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*this->patchNeighbourField()
    );
    Field<Type>::operator=
    (
        intact_*iF + (1.0 - intact_)*intVals
    );
}


template<class Type>
void Foam::burstZeroGradientFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntry(os, "burstP", pBurst_);
    writeEntry(os, "cyclicPatch", cyclicPatchName_);

    if (useImpulse_)
    {
        writeEntry(os, "impulse", impulse_);
        writeEntry(os, "burstImpulse", burstImpulse_);
    }

    writeEntry(os, "partialBurst", partialBurst_);
    writeEntry(os, "intact", intact_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
