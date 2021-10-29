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

#include "temperatureCoupledImmersedMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "thermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

tmp<scalarField> temperatureCoupledImmersedMixedFvPatchScalarField::kappa
(
    const volScalarField::Internal& Tp,
    const label patchi
) const
{
    const fvMesh& mesh = Tp.mesh();
    const word& phase(Tp.group());
    const word thermoName
    (
        IOobject::groupName(basicThermo::dictName, phase)
    );

    if (mesh.foundObject<fluidThermo>(thermoName))
    {
        const word ttmName
        (
            IOobject::groupName
            (
                thermophysicalTransportModel::typeName,
                phase
            )
        );

        if (mesh.foundObject<thermophysicalTransportModel>(ttmName))
        {
            const thermophysicalTransportModel& ttm =
                mesh.lookupObject<thermophysicalTransportModel>(ttmName);

            if (patchi >= 0)
            {
                return ttm.kappaEff(patchi);
            }
            else
            {
                return
                    ibm_.interpolateTo(ttm.kappaEff()().primitiveField());
            }
        }
        else
        {
            const fluidThermo& thermo =
                mesh.lookupObject<fluidThermo>(thermoName);
            if (patchi >= 0)
            {
                return thermo.kappa(patchi);
            }
            else
            {
                return
                    ibm_.interpolateTo(thermo.kappa()().primitiveField());
            }
        }
    }
    else if (mesh.foundObject<solidThermo>(thermoName))
    {
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>(thermoName);

        if (!thermo.isotropic() && patchi >= 0)
        {
            const symmTensorField kappa(thermo.KappaLocal(patchi));
            const vectorField n(this->patch().nf());

            return n & kappa & n;
        }
        else
        {
            if (patchi >= 0)
            {
                return thermo.kappa(patchi);
            }
            else
            {
                return
                    ibm_.interpolateTo(thermo.kappa()().primitiveField());
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find a fluidThermo or solidThermo instance"
            << exit(FatalError);

        return scalarField::null();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

temperatureCoupledImmersedMixedFvPatchScalarField::
temperatureCoupledImmersedMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    name_
    (
        p.boundaryMesh().mesh().dbDir()
    ),
    mapper_
    (
        p.boundaryMesh().mesh().time().lookupObject<immersedMeshMapper>
        (
            IOobject::groupName("immersedMeshMapper", name_)
        )
    ),
    ibm_(mapper_.immersedObject()),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


temperatureCoupledImmersedMixedFvPatchScalarField::
temperatureCoupledImmersedMixedFvPatchScalarField
(
    const temperatureCoupledImmersedMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    name_
    (
        p.boundaryMesh().mesh().dbDir()
    ),
    mapper_
    (
        p.boundaryMesh().mesh().time().lookupObject<immersedMeshMapper>
        (
            IOobject::groupName("immersedMeshMapper", name_)
        )
    ),
    ibm_(mapper_.immersedObject()),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_)
{}


temperatureCoupledImmersedMixedFvPatchScalarField::
temperatureCoupledImmersedMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    name_
    (
        p.boundaryMesh().mesh().dbDir()
    ),
    mapper_
    (
        p.boundaryMesh().mesh().time().lookupObject<immersedMeshMapper>
        (
            IOobject::groupName("immersedMeshMapper", name_)
        )
    ),
    ibm_(mapper_.immersedObject()),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.lookupOrDefault<word>("qrNbr", "none")),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0)
{
//     if (!isA<mappedPatchBase>(this->patch().patch()))
//     {
//         FatalErrorInFunction
//             << "' not type '" << mappedPatchBase::typeName << "'"
//             << "\n    for patch " << p.name()
//             << " of field " << internalField().name()
//             << " in file " << internalField().objectPath()
//             << exit(FatalError);
//     }

    if (dict.found("thicknessLayers"))
    {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
}


temperatureCoupledImmersedMixedFvPatchScalarField::
temperatureCoupledImmersedMixedFvPatchScalarField
(
    const temperatureCoupledImmersedMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    name_(psf.name_),
    mapper_(psf.mapper_),
    ibm_(psf.ibm_),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void temperatureCoupledImmersedMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!ibm_.pMesh().foundObject<volScalarField>(TnbrName_))
    {
        valueFraction() = Zero;
        refValue() = patchInternalField();
        refGrad() = Zero;
        mixedFvPatchScalarField::updateCoeffs();
        return;
    }
    label patchi = this->patch().index();

    const volScalarField& volNbrTp =
        ibm_.pMesh().lookupObject<volScalarField>(TnbrName_);
    scalarField nbrTp
    (
        mapper_.mapObjectToBoundary(ibm_.patchInternalField(volNbrTp)())
    );

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    if (contactRes_ == 0.0)
    {
        KDeltaNbr = kappa(volNbrTp)*ibm_.deltaCoeffs();
    }
    else
    {
        KDeltaNbr.setSize(nbrTp.size(), contactRes_);
    }
    KDeltaNbr = mapper_.mapObjectToBoundary(KDeltaNbr);
    scalarField K(kappa(this->internalField(), patchi));
    scalarField KDelta(K*patch().deltaCoeffs());

    scalarField qr(this->size(), 0.0);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    scalarField qrNbr(this->size(), 0.0);
    if (qrNbrName_ != "none")
    {
        const volScalarField& volQrNbr
        (
            ibm_.pMesh().lookupObject<volScalarField>(qrNbrName_)
        );
        qrNbr = mapper_.mapObjectToBoundary(ibm_.interpolateTo(volQrNbr)());
    }

    valueFraction() = KDeltaNbr/stabilise((KDeltaNbr + KDelta), small);
    refValue() = nbrTp;

    refGrad() = (qr + qrNbr)/max(K, small);

    mixedFvPatchScalarField::updateCoeffs();


//     if (debug)
//     {
//         scalar Q = gSum(K*patch().magSf()*snGrad());
//
//         Info<< patch().boundaryMesh().mesh().name() << ':'
//             << patch().name() << ':'
//             << this->internalField().name() << " <- "
//             << nbrMesh.name() << ':'
//             << nbrPatch.name() << ':'
//             << this->internalField().name() << " :"
//             << " heat transfer rate:" << Q
//             << " walltemperature "
//             << " min:" << gMin(*this)
//             << " max:" << gMax(*this)
//             << " avg:" << gAverage(*this)
//             << endl;
//     }
}


void temperatureCoupledImmersedMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    writeEntry(os, "Tnbr", TnbrName_);
    writeEntry(os, "qrNbr", qrNbrName_);
    writeEntry(os, "qr", qrName_);
    writeEntry(os, "thicknessLayers", thicknessLayers_);
    writeEntry(os, "kappaLayers", kappaLayers_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    temperatureCoupledImmersedMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
