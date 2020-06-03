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

#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "blastFixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastFixedEnergyFvPatchScalarField::
blastFixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::blastFixedEnergyFvPatchScalarField::
blastFixedEnergyFvPatchScalarField
(
    const blastFixedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::blastFixedEnergyFvPatchScalarField::
blastFixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict)
{}


Foam::blastFixedEnergyFvPatchScalarField::
blastFixedEnergyFvPatchScalarField
(
    const blastFixedEnergyFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf)
{}


Foam::blastFixedEnergyFvPatchScalarField::
blastFixedEnergyFvPatchScalarField
(
    const blastFixedEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blastFixedEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    fluidThermoModel& thermo =
        db().lookupObjectRef<fluidThermoModel>
        (
            IOobject::groupName("basicThermo", internalField().group())
        );
    label patchID = patch().index();

    const scalarField& rhow = thermo.rho().boundaryField()[patchID];
    const scalarField& ew = thermo.e().boundaryField()[patchID];
    fvPatchScalarField& Tw =
        const_cast<fvPatchScalarField&>(thermo.T().boundaryFieldRef()[patchID]);

    Tw.evaluate();
    operator==(thermo.e(rhow, ew, Tw, patchID));

    fixedValueFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        blastFixedEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
