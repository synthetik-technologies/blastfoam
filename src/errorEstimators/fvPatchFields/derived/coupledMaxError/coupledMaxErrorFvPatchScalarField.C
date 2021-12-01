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
#include "errorEstimator.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const coupledMaxErrorFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{}


Foam::coupledMaxErrorFvPatchScalarField::coupledMaxErrorFvPatchScalarField
(
    const coupledMaxErrorFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledMaxErrorFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    fixedValueFvPatchField<scalar>::updateCoeffs();
    Field<scalar>::operator=(this->patchInternalField());

    patch().boundaryMesh().mesh().lookupObjectRef<errorEstimator>
    (
        errorEstimator::typeName
    ).update();


    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(this->patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());

    scalarField::operator=(this->patchInternalField());

    if (nbrMesh.foundObject<errorEstimator>(errorEstimator::typeName))
    {
        errorEstimator& nbrErrorEst =
            nbrMesh.lookupObjectRef<errorEstimator>
            (
                errorEstimator::typeName
            );

        nbrErrorEst.update();

        const label samplePatchi = mpp.samplePolyPatch().index();
        scalarField nbrError
        (
            nbrErrorEst.error().boundaryField()
            [
                samplePatchi
            ].patchInternalField()
        );

        mpp.distribute(nbrError);

        scalarField::operator=(max(nbrError, *this));

        volScalarField::Internal& iF
        (
            const_cast<volScalarField::Internal&>(this->internalField())
        );
        forAll(*this, facei)
        {
            iF[patch().faceCells()[facei]] = this->operator[](facei);
        }
    }
}


void Foam::coupledMaxErrorFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
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
