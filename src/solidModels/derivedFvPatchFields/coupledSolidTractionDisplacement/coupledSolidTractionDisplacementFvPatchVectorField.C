/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
25-06-2021 Synthetik Applied Technologies: |    Added coupledSolidTraction
-------------------------------------------------------------------------------
License
    This file is a derivatived work of OpenFOAM.

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

#include "coupledSolidTractionDisplacementFvPatchVectorField.H"
#include "uniformDimensionedFields.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledSolidTractionDisplacementFvPatchVectorField::
coupledSolidTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    mpp_(p),
    pName_("p"),
    pRef_(0.0)
{}


Foam::coupledSolidTractionDisplacementFvPatchVectorField::
coupledSolidTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    mpp_(p),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_(dict.lookupOrDefault("pRef", 0.0))
{}


Foam::coupledSolidTractionDisplacementFvPatchVectorField::
coupledSolidTractionDisplacementFvPatchVectorField
(
    const coupledSolidTractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, p, iF, mapper),
    mpp_(p),
    pName_(tdpvf.pName_),
    pRef_(tdpvf.pRef_)
{}


Foam::coupledSolidTractionDisplacementFvPatchVectorField::
coupledSolidTractionDisplacementFvPatchVectorField
(
    const coupledSolidTractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(tdpvf, iF),
    mpp_(tdpvf.mpp_),
    pName_(tdpvf.pName_),
    pRef_(tdpvf.pRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledSolidTractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void Foam::coupledSolidTractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}

void Foam::coupledSolidTractionDisplacementFvPatchVectorField::updateCoeffs()
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
    const polyMesh& nbrMesh = mpp_.sampleMesh();
    const label samplePatchi = mpp_.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    //- Lookup viscous stress and pressure fields
    symmTensorField nbrDevRhoReff
    (
        nbrPatch.lookupPatchField<volSymmTensorField, symmTensor>("devRhoReff")
    );
    scalarField nbrP
    (
        nbrPatch.lookupPatchField<volScalarField, scalar>(pName_)
    );

    mpp_.distribute(nbrDevRhoReff);
    mpp_.distribute(nbrP);

    this->pressure() = nbrP - pRef_;
    this->traction() = nbrDevRhoReff & patch().nf();

    tractionDisplacementFvPatchVectorField::updateCoeffs();
}


void Foam::coupledSolidTractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);
    writeEntry(os, "pName", pName_);
    writeEntry(os, "pRef", pRef_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        coupledSolidTractionDisplacementFvPatchVectorField
    );
}


// ************************************************************************* //
