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

#include "blastMixedEnergyFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermoModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastMixedEnergyFvPatchScalarField::
blastMixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF)
{
    valueFraction() = 0.0;
    refValue() = 0.0;
    refGrad() = 0.0;
}


Foam::blastMixedEnergyFvPatchScalarField::
blastMixedEnergyFvPatchScalarField
(
    const blastMixedEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::blastMixedEnergyFvPatchScalarField::
blastMixedEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict)
{}


Foam::blastMixedEnergyFvPatchScalarField::
blastMixedEnergyFvPatchScalarField
(
    const blastMixedEnergyFvPatchScalarField& tppsf
)
:
    mixedFvPatchScalarField(tppsf)
{}


Foam::blastMixedEnergyFvPatchScalarField::
blastMixedEnergyFvPatchScalarField
(
    const blastMixedEnergyFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blastMixedEnergyFvPatchScalarField::updateCoeffs()
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
    mixedFvPatchScalarField& Tw = refCast<mixedFvPatchScalarField>
    (
        const_cast<fvPatchScalarField&>(thermo.T().boundaryField()[patchID])
    );

    Tw.evaluate();

    valueFraction() = Tw.valueFraction();
    refValue() = thermo.e(rhow, ew, Tw.refValue(), patchID);
    refGrad() =
        thermo.Cv(rhow, ew, Tw, patchID)*Tw.refGrad()
      + patch().deltaCoeffs()*
        (
            thermo.e(rhow, ew, Tw, patchID)
          - thermo.e(rhow, ew, Tw, patch().faceCells())
        );

    mixedFvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        blastMixedEnergyFvPatchScalarField
    );
}

// ************************************************************************* //
