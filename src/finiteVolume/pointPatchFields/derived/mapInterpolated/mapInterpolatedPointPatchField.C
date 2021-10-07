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

#include "mapInterpolatedPointPatchField.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mapInterpolatedPointPatchField<Type>::mapInterpolatedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    mpp_
    (
        dynamicCast<const fvMesh>
        (
            p.boundaryMesh().mesh().thisDb()
        ).boundary()[p.index()]
    ),
    volName_(word(iF.name()).replaceAll("point", word::null))
{}


template<class Type>
Foam::mapInterpolatedPointPatchField<Type>::mapInterpolatedPointPatchField
(
    const mapInterpolatedPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(p, iF),
    mpp_
    (
        dynamicCast<const fvMesh>
        (
            p.boundaryMesh().mesh().thisDb()
        ).boundary()[p.index()]
    ),
    volName_(ptf.volName_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        pointPatchField<Type>::operator=(this->patchInternalField());
    }
    mapper(*this, ptf);
}


template<class Type>
Foam::mapInterpolatedPointPatchField<Type>::mapInterpolatedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    mpp_
    (
        dynamicCast<const fvMesh>
        (
            p.boundaryMesh().mesh().thisDb()
        ).boundary()[p.index()]
    ),
    volName_
    (
        dict.lookupOrDefault<word>
        (
            "volName",
            word(iF.name()).replaceAll("point", word::null)
        )
    )
{}


template<class Type>
Foam::mapInterpolatedPointPatchField<Type>::mapInterpolatedPointPatchField
(
    const mapInterpolatedPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    mpp_(ptf.mpp_),
    volName_(ptf.volName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mapInterpolatedPointPatchField<Type>::updateCoeffs()
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
    const polyPatch& pPatch =
        dynamicCast<const polyMesh>
        (
            this->patch().boundaryMesh().mesh().thisDb()
        ).boundaryMesh()[this->patch().index()];
    const polyMesh& nbrMesh = mpp_.sampleMesh();
    const polyPatch& samplePatch = mpp_.samplePolyPatch();
    const label samplePatchi = samplePatch.index();

    const GeometricField<Type, fvPatchField, volMesh>& nbr =
        nbrMesh.lookupObject<GeometricField<Type, fvPatchField, volMesh>>(volName_);
    Field<Type> pNbr(nbr.boundaryField()[samplePatchi]);

    mpp_.distribute(pNbr);
    primitivePatchInterpolation ppI(pPatch);

    Field<Type>::operator=(ppI.faceToPointInterpolate(pNbr));

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mapInterpolatedPointPatchField<Type>::write(Ostream& os) const
{
    pointPatchField<Type>::write(os);
    writeEntry(os, "volName", volName_);
}

// ************************************************************************* //
