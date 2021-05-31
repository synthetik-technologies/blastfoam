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

#include "coupledMotionPointPatchField.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "mappedPatchBase.H"
#include "mappedMovingPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::coupledMotionPointPatchField<Type>::coupledMotionPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    mpp_(p),
    DnbrName_("pointD"),
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
Foam::coupledMotionPointPatchField<Type>::coupledMotionPointPatchField
(
    const coupledMotionPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(p, iF),
    mpp_(p),
    DnbrName_(ptf.DnbrName_),
    cmpt_(ptf.cmpt_),
    velocity_(ptf.velocity_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        pointPatchField<Type>::operator=(this->patchInternalField());
    }
    mapper(*this, ptf);
}


template<class Type>
Foam::coupledMotionPointPatchField<Type>::coupledMotionPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    mpp_(p),
    DnbrName_(dict.lookupOrDefault<word>("DnbrName", "pointD")),
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
Foam::coupledMotionPointPatchField<Type>::coupledMotionPointPatchField
(
    const coupledMotionPointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    mpp_(ptf.mpp_),
    DnbrName_(ptf.DnbrName_),
    cmpt_(ptf.cmpt_),
    velocity_(ptf.velocity_)
{}


template<class Type>
Foam::coupledMotionPointPatchField<Type>::coupledMotionPointPatchField
(
    const coupledMotionPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    mpp_(ptf.mpp_),
    DnbrName_(ptf.DnbrName_),
    cmpt_(ptf.cmpt_),
    velocity_(ptf.velocity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::coupledMotionPointPatchField<Type>::updateCoeffs()
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

    const pointVectorField& nbrD =
        nbrFvMesh.lookupObject<pointVectorField>(DnbrName_);
    vectorField Dp
    (
        nbrD.boundaryField()[samplePatchi].patchInternalField()
    );

    mpp_.distribute(Dp);

    if (velocity_)
    {
        vectorField DpOld
        (
            nbrD.oldTime().boundaryField()[samplePatchi].patchInternalField()
        );
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

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::coupledMotionPointPatchField<Type>::write(Ostream& os) const
{
    pointPatchField<Type>::write(os);
    writeEntry(os, "DnbrName", DnbrName_);
}

// ************************************************************************* //
