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

#include "globalMappedPointPatchField.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "coupledGlobalPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(iF.name())
{}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const globalMappedPointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(p, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(ptf.nbrName_)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        pointPatchField<Type>::operator=(this->patchInternalField());
    }
    mapper(*this, ptf);
}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                p.boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(dict.lookup<word>("nbrName"))
{}


template<class Type>
Foam::globalMappedPointPatchField<Type>::globalMappedPointPatchField
(
    const globalMappedPointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New
        (
            dynamicCast<const polyMesh>
            (
                ptf.patch().boundaryMesh().mesh().thisDb()
            )
        )
    ),
    nbrName_(ptf.nbrName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::globalMappedPointPatchField<Type>::updateCoeffs()
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
    const coupledGlobalPolyPatch& cgpp = globalBoundary_(this->patch());
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const label samplePatchi = samplePatch.patch().index();

    const pointPatchField<Type>& nbr =
        nbrMesh.lookupObject<GeometricField<Type, pointPatchField, pointMesh>>
        (
            nbrName_
        ).boundaryField()[samplePatchi];

    Field<Type>::operator=
    (
        samplePatch.pointInterpolate(nbr.patchInternalField())
    );
    fixedValuePointPatchField<Type>::updateCoeffs();

    // Restore tag
    UPstream::msgType() = oldTag;
}


template<class Type>
void Foam::globalMappedPointPatchField<Type>::write(Ostream& os) const
{
    fixedValuePointPatchField<Type>::write(os);
    writeEntry(os, "nbrName", nbrName_);
}

// ************************************************************************* //
