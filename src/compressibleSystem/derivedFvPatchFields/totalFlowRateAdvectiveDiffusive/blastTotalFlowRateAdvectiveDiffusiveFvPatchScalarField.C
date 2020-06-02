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

#include "blastTotalFlowRateAdvectiveDiffusiveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOobjectList.H"
#include "blastCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    rhoName_("none"),
    massFluxFraction_(1.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "none")),
    massFluxFraction_(dict.lookupOrDefault<scalar>("massFluxFraction", 1.0))
{

    refValue() = 1.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            Field<scalar>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}

Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    massFluxFraction_(ptf.massFluxFraction_)
{}


Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}

Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
totalFlowRateAdvectiveDiffusiveFvPatchScalarField
(
    const totalFlowRateAdvectiveDiffusiveFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    massFluxFraction_(tppsf.massFluxFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    m(*this, *this);
}


void Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);
}


void
Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const label patchi = patch().index();

    const blastCompressibleTurbulenceModel& turbModel =
        db().lookupObject<blastCompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    const scalarField alphap(turbModel.alphaEff(patchi));

    refValue() = massFluxFraction_;
    refGrad() = 0.0;

    valueFraction() =
        1.0
        /
        (
            1.0 +
            alphap*patch().deltaCoeffs()*patch().magSf()/max(mag(phip), small)
        );

    mixedFvPatchField<scalar>::updateCoeffs();

    if (debug)
    {
        scalar phi = gSum(-phip*(*this));

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " mass flux[Kg/s]:" << phi
            << endl;
    }
}


void Foam::blast::totalFlowRateAdvectiveDiffusiveFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry(os, "phi", phiName_);
    writeEntry(os, "rho", rhoName_);
    writeEntry(os, "massFluxFraction", massFluxFraction_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blast
{
    makePatchTypeField
    (
        fvPatchScalarField,
        totalFlowRateAdvectiveDiffusiveFvPatchScalarField
    );

}
}

// ************************************************************************* //
