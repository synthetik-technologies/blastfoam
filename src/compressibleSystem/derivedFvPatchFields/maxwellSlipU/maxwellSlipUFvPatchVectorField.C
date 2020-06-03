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

#include "maxwellSlipUFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "fluidThermoModel.H"
#include "thermodynamicConstants.H"
#include "blastCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    accommodationCoeff_(1.0),
    Uwall_(p.size(), vector(0.0, 0.0, 0.0)),
    thermalCreep_(true),
    curvature_(true)
{}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& mspvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFixedValueSlipFvPatchVectorField(mspvf, p, iF, mapper),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFixedValueSlipFvPatchVectorField(p, iF),
    accommodationCoeff_(readScalar(dict.lookup("accommodationCoeff"))),
    Uwall_("Uwall", dict, p.size()),
    thermalCreep_(dict.lookupOrDefault("thermalCreep", true)),
    curvature_(dict.lookupOrDefault("curvature", true))
{
    if
    (
        mag(accommodationCoeff_) < small
     || mag(accommodationCoeff_) > 2.0
    )
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "unphysical accommodationCoeff_ specified"
            << "(0 < accommodationCoeff_ <= 1)" << endl
            << exit(FatalIOError);
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );

        if (dict.found("refValue") && dict.found("valueFraction"))
        {
            this->refValue() = vectorField("refValue", dict, p.size());
            this->valueFraction() =
                scalarField("valueFraction", dict, p.size());
        }
        else
        {
            this->refValue() = *this;
            this->valueFraction() = scalar(1);
        }
    }
}


Foam::maxwellSlipUFvPatchVectorField::maxwellSlipUFvPatchVectorField
(
    const maxwellSlipUFvPatchVectorField& mspvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFixedValueSlipFvPatchVectorField(mspvf, iF),
    accommodationCoeff_(mspvf.accommodationCoeff_),
    Uwall_(mspvf.Uwall_),
    thermalCreep_(mspvf.thermalCreep_),
    curvature_(mspvf.curvature_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::maxwellSlipUFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fluidThermoModel& thermo =
        db().lookupObject<fluidThermoModel>
        (
            IOobject::groupName("basicThermo", internalField().group())
        );
    const label patchi = patch().index();
    const scalarField& pmu = thermo.mu(patchi);
    const scalarField& prho = thermo.rho().boundaryField()[patchi];
    const volScalarField& vsfT = thermo.T();
    const fvPatchScalarField& pT = vsfT.boundaryField()[patchi];

    Field<scalar> C1
    (
        sqrt
        (
            Foam::constant::thermodynamic::RR/thermo.W(patchi)*pT
           *constant::mathematical::piByTwo)
      * (2.0 - accommodationCoeff_)/accommodationCoeff_
    );

    Field<scalar> pnu(pmu/prho);
    valueFraction() = (1.0/(1.0 + patch().deltaCoeffs()*C1*pnu));

    refValue() = Uwall_;

    if (thermalCreep_)
    {

        Field<vector> gradpT(fvc::grad(vsfT)().boundaryField()[patchi]);
        vectorField n(patch().nf());

        refValue() -= 3.0*pnu/(4.0*pT)*transform(I - n*n, gradpT);
    }

    if (curvature_)
    {
        const volVectorField& vsfU =
            db().lookupObject<volVectorField>(internalField().name());
        const blastCompressibleTurbulenceModel& turbulence =
            db().lookupObject<blastCompressibleTurbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    internalField().group()
                )
            );
        tensorField ptauMC
        (
            turbulence.muEff(patchi)
           *dev2
            (
                Foam::T(fvc::grad(vsfU))().boundaryField()[patchi]
            )
        );
        vectorField n(patch().nf());

        refValue() -= C1/prho*transform(I - n*n, (n & ptauMC));
    }

    mixedFixedValueSlipFvPatchVectorField::updateCoeffs();
}


void Foam::maxwellSlipUFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    writeEntry(os, "accommodationCoeff", accommodationCoeff_);
    writeEntry(os, "Uwall", Uwall_);
    writeEntry(os, "thermalCreep", thermalCreep_);
    writeEntry(os, "curvature", curvature_);

    writeEntry(os, "refValue", refValue());
    writeEntry(os, "valueFraction", valueFraction());

    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        maxwellSlipUFvPatchVectorField
    );
}

// ************************************************************************* //
