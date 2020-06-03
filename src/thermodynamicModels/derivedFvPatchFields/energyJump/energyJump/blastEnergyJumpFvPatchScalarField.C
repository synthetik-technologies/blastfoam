/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "addToRunTimeSelectionTable.H"
#include "blastEnergyJumpFvPatchScalarField.H"
#include "fixedJumpFvPatchFields.H"
#include "fluidThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastEnergyJumpFvPatchScalarField::blastEnergyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{}


Foam::blastEnergyJumpFvPatchScalarField::blastEnergyJumpFvPatchScalarField
(
    const blastEnergyJumpFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<scalar>(ptf, p, iF, mapper)
{}


Foam::blastEnergyJumpFvPatchScalarField::blastEnergyJumpFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<scalar>(p, iF)
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


Foam::blastEnergyJumpFvPatchScalarField::blastEnergyJumpFvPatchScalarField
(
    const blastEnergyJumpFvPatchScalarField& ptf
)
:
    fixedJumpFvPatchField<scalar>(ptf)
{}


Foam::blastEnergyJumpFvPatchScalarField::blastEnergyJumpFvPatchScalarField
(
    const blastEnergyJumpFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpFvPatchField<scalar>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blastEnergyJumpFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        fluidThermoModel& thermo =
            db().lookupObjectRef<fluidThermoModel>
            (
                IOobject::groupName("basicThermo", internalField().group())
            );
        label patchID = patch().index();

        const scalarField& rho = thermo.rho().boundaryField()[patchID];
        const scalarField& e = thermo.e().boundaryField()[patchID];
        fvPatchScalarField& TbPatch = thermo.T().boundaryFieldRef()[patchID];

        fixedJumpFvPatchScalarField& Tbp =
            refCast<fixedJumpFvPatchScalarField>(TbPatch);

        // force update of jump
        Tbp.updateCoeffs();

        const labelUList& faceCells = this->patch().faceCells();

        jump_ = thermo.e(rho, e, Tbp.jump(), faceCells);
    }

    fixedJumpFvPatchField<scalar>::updateCoeffs();
}


void Foam::blastEnergyJumpFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpFvPatchField<scalar>::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       blastEnergyJumpFvPatchScalarField
   );
}

// ************************************************************************* //
