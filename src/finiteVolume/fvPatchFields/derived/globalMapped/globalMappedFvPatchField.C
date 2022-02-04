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

#include "globalMappedFvPatchField.H"
#include "coupledGlobalPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::globalMappedFvPatchField<Type>::globalMappedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    globalBoundary_(globalPolyBoundaryMesh::New(p.boundaryMesh().mesh())),
    nbrName_(iF.name())
{}


template<class Type>
Foam::globalMappedFvPatchField<Type>::globalMappedFvPatchField
(
    const globalMappedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    globalBoundary_(globalPolyBoundaryMesh::New(p.boundaryMesh().mesh())),
    nbrName_(iF.name())
{}


template<class Type>
Foam::globalMappedFvPatchField<Type>::globalMappedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    globalBoundary_(globalPolyBoundaryMesh::New(p.boundaryMesh().mesh())),
    nbrName_(dict.lookup<word>("nbrName"))
{
    fvPatchField<Type>::operator==(this->patchInternalField());
}


template<class Type>
Foam::globalMappedFvPatchField<Type>::globalMappedFvPatchField
(
    const globalMappedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New(ptf.patch().boundaryMesh().mesh())
    ),
    nbrName_(iF.name())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::globalMappedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
}


template<class Type>
void Foam::globalMappedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::globalMappedFvPatchField<Type>::patchNeighbourField() const
{
    const coupledGlobalPolyPatch& cgpp =
        globalBoundary_(this->patch().patch());
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const label samplePatchi = samplePatch.patch().index();

    const fvPatchField<Type>& nbr =
        nbrMesh.lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            nbrName_
        ).boundaryField()[samplePatchi];
    return cgpp.faceInterpolate(nbr.patchInternalField());

}


template<class Type>
void Foam::globalMappedFvPatchField<Type>::updateCoeffs()
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
    const coupledGlobalPolyPatch& cgpp =
        globalBoundary_(this->patch().patch());
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const label samplePatchi = samplePatch.patch().index();

    const fvPatchField<Type>& nbr =
        nbrMesh.lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            nbrName_
        ).boundaryField()[samplePatchi];


    Field<Type>::operator==(samplePatch.pointInterpolate(nbr));
    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::globalMappedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "nbrName", nbrName_);
}


// ************************************************************************* //
