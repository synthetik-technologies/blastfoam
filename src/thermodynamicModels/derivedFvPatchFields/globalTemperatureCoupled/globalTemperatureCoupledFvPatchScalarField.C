/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "globalTemperatureCoupledFvPatchScalarField.H"
#include "solidThermo.H"
#include "fluidThermo.H"
#include "thermophysicalTransportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "coupledGlobalPolyPatch.H"
#include "thermodynamicConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

Foam::tmp<Foam::scalarField>
Foam::globalTemperatureCoupledFvPatchScalarField::kappa
(
    const fvPatchScalarField& Tp
)
{
    const fvMesh& mesh = Tp.patch().boundaryMesh().mesh();
    const label patchi = Tp.patch().index();

    const word& phase(Tp.internalField().group());

    const word thermoName
    (
        IOobject::groupName(basicThermo::dictName, phase)
    );

    if (mesh.foundObject<fluidThermo>(thermoName))
    {
        static word ttmName
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

            return ttm.kappaEff(patchi);
        }
        else
        {
            const fluidThermo& thermo =
                mesh.lookupObject<fluidThermo>(thermoName);

            return thermo.kappa(patchi);
        }
    }
    else if (mesh.foundObject<solidThermo>(thermoName))
    {
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>(thermoName);

        if (!thermo.isotropic())
        {
            const symmTensorField kappa(thermo.KappaLocal(patchi));
            const vectorField n(Tp.patch().nf());

            return n & kappa & n;
        }
        else
        {
            return thermo.kappa(patchi);
        }
    }
//     else if
//     (
//         mesh.foundObject<volScalarField>
//         (
//             IOobject::groupName("kappa", phase)
//         )
//     )
    {
        return Tp.patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("kappa", phase)
        );
    }
//     else
//     {
//         FatalErrorInFunction
//             << "Cannot find a fluidThermo or solidThermo instance"
//             << exit(FatalError);
//
//         return scalarField::null();
//     }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

globalTemperatureCoupledFvPatchScalarField::
globalTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    globalBoundary_(globalPolyBoundaryMesh::New(p.boundaryMesh().mesh())),
    TnbrName_("T"),
    hNbrName_("none"),
    hName_("none"),
    qrNbrName_("none"),
    qrName_("none"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),
    unmappedT_(constant::thermodynamic::Tstd)
{
    this->refValue() = unmappedT_;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


globalTemperatureCoupledFvPatchScalarField::
globalTemperatureCoupledFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    globalBoundary_(globalPolyBoundaryMesh::New(p.boundaryMesh().mesh())),
    TnbrName_(dict.lookupOrDefault<word>("TNbr", "T")),
    hNbrName_(dict.lookupOrDefault<word>("hNbr", "none")),
    hName_(dict.lookupOrDefault<word>("h", "none")),
    qrNbrName_(dict.lookupOrDefault<word>("qrNbr", "none")),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),
    unmappedT_
    (
        dict.lookupOrDefault("unmappedT", constant::thermodynamic::Tstd)
    )
{
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


globalTemperatureCoupledFvPatchScalarField::
globalTemperatureCoupledFvPatchScalarField
(
    const globalTemperatureCoupledFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    globalBoundary_(globalPolyBoundaryMesh::New(p.boundaryMesh().mesh())),
    TnbrName_(psf.TnbrName_),
    hNbrName_(psf.hNbrName_),
    hName_(psf.hName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_),
    unmappedT_(psf.unmappedT_)
{}


globalTemperatureCoupledFvPatchScalarField::
globalTemperatureCoupledFvPatchScalarField
(
    const globalTemperatureCoupledFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    globalBoundary_
    (
        globalPolyBoundaryMesh::New(psf.patch().boundaryMesh().mesh())
    ),
    TnbrName_(psf.TnbrName_),
    hNbrName_(psf.hNbrName_),
    hName_(psf.hName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_),
    unmappedT_(psf.unmappedT_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void globalTemperatureCoupledFvPatchScalarField::updateCoeffs()
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
    const coupledGlobalPolyPatch& cgpp =
        globalBoundary_(this->patch().patch());
    const polyMesh& nbrMesh = cgpp.sampleMesh();
    const coupledGlobalPolyPatch& samplePatch = cgpp.samplePatch();
    const label samplePatchi = samplePatch.patch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    if (!nbrMesh.foundObject<volScalarField>(TnbrName_))
    {
        if (debug)
        {
            WarningInFunction
                << TnbrName_ << " was not found in " << nbrMesh.name()
                << endl;
        }
        mixedFvPatchScalarField::updateCoeffs();

        // Restore tag
        UPstream::msgType() = oldTag;
        return;
    }

    scalarField TcOwn(patchInternalField());
    scalarField& Tp = *this;

    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);

    // Swap to obtain full local values of neighbour internal field
    scalarField TcNbr
    (
        samplePatch.faceInterpolate(nbrTp.patchInternalField())
    );

    //- Difference in temperature
    const scalarField deltaT(TcNbr - TcOwn);

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    if (contactRes_ == 0.0)
    {
        KDeltaNbr =
            samplePatch.faceInterpolate
            (
                kappa(nbrTp)*nbrPatch.deltaCoeffs()
            );
    }
    else
    {
        KDeltaNbr.setSize(this->size(), contactRes_);
    }

    scalarField KDelta(kappa(*this)*patch().deltaCoeffs());

    scalarField q(Tp.size(), 0.0);
    if (hName_ != "none")
    {
        q +=
            patch().lookupPatchField<volScalarField, scalar>(hName_)
           *deltaT;
    }

    if (hNbrName_ != "none")
    {
        q +=
            samplePatch.faceInterpolate
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (
                    hNbrName_
                )
            )*deltaT;
    }

    if (qrName_ != "none")
    {
        q += patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    if (qrNbrName_ != "none")
    {
        q +=
            samplePatch.faceInterpolate
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (
                    qrNbrName_
                )
            );
    }

    valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
    refValue() = TcNbr;
    refGrad() = q/kappa(*this);

    if (cgpp.hasUnmappedFaces())
    {
        cgpp.setUnmappedFace(valueFraction(), 1.0);
        cgpp.setUnmappedFace(refValue(), unmappedT_);
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(kappa(*this)*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << " walltemperature "
            << " min:" << gMin(Tp)
            << " max:" << gMax(Tp)
            << " avg:" << gAverage(Tp)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void globalTemperatureCoupledFvPatchScalarField::write
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
    writeEntry(os, "unmappedT", unmappedT_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    globalTemperatureCoupledFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
