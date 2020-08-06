/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "FSIForceDisplacementFvPatchVectorField.H"
#include "uniformDimensionedFields.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSIForceDisplacementFvPatchVectorField::
FSIForceDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    pName_("p"),
    pRef_(this->size(), 0.0)
{}


Foam::FSIForceDisplacementFvPatchVectorField::
FSIForceDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_("pRef", dict, p.size())
{}


Foam::FSIForceDisplacementFvPatchVectorField::
FSIForceDisplacementFvPatchVectorField
(
    const FSIForceDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, p, iF, mapper),
    pName_(tdpvf.pName_),
    pRef_(mapper(tdpvf.pRef_))
{}


Foam::FSIForceDisplacementFvPatchVectorField::
FSIForceDisplacementFvPatchVectorField
(
    const FSIForceDisplacementFvPatchVectorField& tdpvf
)
:
    tractionDisplacementFvPatchVectorField(tdpvf),
    pName_(tdpvf.pName_),
    pRef_(tdpvf.pRef_)
{}


Foam::FSIForceDisplacementFvPatchVectorField::
FSIForceDisplacementFvPatchVectorField
(
    const FSIForceDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, iF),
    pName_(tdpvf.pName_),
    pRef_(tdpvf.pRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FSIForceDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
    m(pRef_, pRef_);
}


void Foam::FSIForceDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);

    const FSIForceDisplacementFvPatchVectorField& dmptf =
        refCast<const FSIForceDisplacementFvPatchVectorField>(ptf);

    pRef_.rmap(dmptf.pRef_, addr);
}

void Foam::FSIForceDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    //- Lookup viscous stress and pressure fields
    vectorField normal(nbrPatch.Sf()/nbrPatch.magSf());
    symmTensorField nbrDevRhoReff
    (
        nbrPatch.lookupPatchField<volSymmTensorField, symmTensor>("devRhoReff")
    );
    scalarField nbrP
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>(pName_)
    );

    mpp.distribute(normal);
    mpp.distribute(nbrDevRhoReff);
    mpp.distribute(nbrP);

    this->pressure() = nbrP - pRef_;
    this->traction() = normal & nbrDevRhoReff;

    tractionDisplacementFvPatchVectorField::updateCoeffs();
}


void Foam::FSIForceDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "pName", pName_);
    writeEntry(os, "pRef", pRef_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        FSIForceDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
