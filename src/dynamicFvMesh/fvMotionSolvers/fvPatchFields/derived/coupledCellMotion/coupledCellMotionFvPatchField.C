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

#include "coupledCellMotionFvPatchField.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "mappedMovingPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::coupledCellMotionFvPatchField<Type>::coupledCellMotionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mpp_(p),
    DnbrName_("D"),
    cmpt_(-1),
    velocity_(this->internalField().dimensions() == dimVelocity)
{
    if (pTraits<Type>::nComponents == 1)
    {
        word fieldName = this->internalField().name();
        char cmptName = fieldName.back();
        cmpt_ = cmpt(word(cmptName));
    }
    else if (pTraits<Type>::nComponents != 3)
    {
        FatalErrorInFunction
            << "Only scalar or vector fields are supported, but type " << nl
            << pTraits<Type>::typeName << " was specified." << nl
            << abort(FatalError);
    }
}


template<class Type>
Foam::coupledCellMotionFvPatchField<Type>::coupledCellMotionFvPatchField
(
    const coupledCellMotionFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mpp_(p),
    DnbrName_(ptf.DnbrName_),
    cmpt_(ptf.cmpt_),
    velocity_(ptf.velocity_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }
    mapper(*this, ptf);
}


template<class Type>
Foam::coupledCellMotionFvPatchField<Type>::coupledCellMotionFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    mpp_(p),
    DnbrName_(dict.lookup("DnbrName")),
    cmpt_(-1),
    velocity_(this->internalField().dimensions() == dimVelocity)
{
    if (pTraits<Type>::nComponents == 1)
    {
        word fieldName = this->internalField().name();
        char cmptName = fieldName.back();
        cmpt_ = cmpt(word(cmptName));
    }
    else if (pTraits<Type>::nComponents != 3)
    {
        FatalErrorInFunction
            << "Only scalar or vector fields are supported, but type "
            << pTraits<Type>::typeName << " was specified." << nl
            << abort(FatalError);
    }
}


template<class Type>
Foam::coupledCellMotionFvPatchField<Type>::coupledCellMotionFvPatchField
(
    const coupledCellMotionFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    mpp_(ptf.mpp_),
    DnbrName_(ptf.DnbrName_),
    cmpt_(ptf.cmpt_),
    velocity_(ptf.velocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::coupledCellMotionFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const polyMesh& nbrMesh = mpp_.sampleMesh();
    const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);
    const label samplePatchi = mpp_.samplePolyPatch().index();

    const volVectorField& nbrD =
        nbrFvMesh.lookupObject<volVectorField>(DnbrName_);
    vectorField Dp(nbrD.boundaryField()[samplePatchi]);

    mpp_.distribute(Dp);

    if (velocity_)
    {
        vectorField DpOld(nbrD.oldTime().boundaryField()[samplePatchi]);
        mpp_.distribute(DpOld);
        Field<Type>::operator=
        (
            (getCmpt(Dp) - getCmpt(DpOld))/nbrFvMesh.time().deltaTValue()
        );
    }
    else
    {
        Field<Type>::operator=(getCmpt(Dp));
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::coupledCellMotionFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "value", *this);
    writeEntry(os, "DnbrName", DnbrName_);
}

// ************************************************************************* //
