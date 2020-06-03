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

#include "blastGradientEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastGradientEnergyFvPatchScalarField::
blastGradientEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF)
{}


Foam::blastGradientEnergyFvPatchScalarField::
blastGradientEnergyFvPatchScalarField
(
    const blastGradientEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::blastGradientEnergyFvPatchScalarField::
blastGradientEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict)
{}


Foam::blastGradientEnergyFvPatchScalarField::
blastGradientEnergyFvPatchScalarField
(
    const blastGradientEnergyFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf)
{}


Foam::blastGradientEnergyFvPatchScalarField::
blastGradientEnergyFvPatchScalarField
(
    const blastGradientEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blastGradientEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    basicThermoModel& thermo
    (
        db().lookupObjectRef<basicThermoModel>
        (
            IOobject::groupName("basicThermo", internalField().group())
        )
    );

    label patchID = patch().index();

    const scalarField& rhow = thermo.rho().boundaryField()[patchID];
    const scalarField& ew = thermo.e().boundaryField()[patchID];
    fvPatchScalarField& Tw = thermo.T().boundaryFieldRef()[patchID];

    gradient() = thermo.Cv(rhow, ew, Tw, patchID)*Tw.snGrad()
      + patch().deltaCoeffs()*
        (
            thermo.e(rhow, ew, Tw, patchID)
          - thermo.e(rhow, ew, Tw, patch().faceCells())
        );

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::blastGradientEnergyFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        blastGradientEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
