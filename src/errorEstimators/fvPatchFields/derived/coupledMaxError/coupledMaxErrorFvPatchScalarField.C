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

#include "coupledMaxErrorFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "volFields.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    mpp_(p)
{}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const coupledMaxErrorFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    mpp_(p)
{
    // For unmapped faces set to internal field value (zero-gradient)
    if (notNull(iF) && mapper.hasUnmapped())
    {
        fvPatchField<scalar>::operator=(this->patchInternalField());
    }
    mapper(*this, ptf);
}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    mpp_(p)
{}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const coupledMaxErrorFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    mpp_(ptf.mpp_)
{}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const coupledMaxErrorFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    mpp_(ptf.mpp_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMaxErrorFvPatchScalarField::updateCoeffs()
{
//     if (this->updated())
//     {
//         return;
//     }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const polyMesh& nbrMesh = mpp_.sampleMesh();
    const fvMesh& nbrFvMesh = refCast<const fvMesh>(nbrMesh);
    const label samplePatchi = mpp_.samplePolyPatch().index();
    const fvPatch& nbrPatch = nbrFvMesh.boundary()[samplePatchi];

    scalarField::operator=(this->patchInternalField());

    if (nbrFvMesh.foundObject<volScalarField>("error"))
    {
        const fvPatchScalarField& nbrError
        (
            nbrPatch.lookupPatchField<volScalarField, scalar>("error")
        );
        scalarField errorp(nbrError.patchInternalField());

        mpp_.distribute(errorp);

        scalarField::operator=(max(errorp, *this));

        volScalarField& error
        (
            this->db().lookupObjectRef<volScalarField>(this->internalField().name())
        );
        forAll(*this, facei)
        {
            error[patch().faceCells()[facei]] = this->operator[](facei);
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void Foam::coupledMaxErrorFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        coupledMaxErrorFvPatchScalarField
    );
}

// ************************************************************************* //
