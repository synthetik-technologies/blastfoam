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

template<>
const char* NamedEnum<globalTemperatureCoupledFvPatchScalarField::TRefType, 3>::names[] =
{
    "max",
    "neighbour",
    "kappaByDelta"
};

const NamedEnum<globalTemperatureCoupledFvPatchScalarField::TRefType, 3>
    globalTemperatureCoupledFvPatchScalarField::TRefTypeNames_;

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
        IOobject::groupName(physicalProperties::typeName, phase)
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

            return thermo.kappa().boundaryField()[patchi];
        }
    }
    else if (mesh.foundObject<solidThermo>(thermoName))
    {
        const solidThermo& thermo =
            mesh.lookupObject<solidThermo>(thermoName);

        return thermo.kappa().boundaryField()[patchi];
    }
    else if
    (
        mesh.foundObject<volScalarField>
        (
            IOobject::groupName("kappa", phase)
        )
    )
    {
        return Tp.patch().lookupPatchField<volScalarField, scalar>
        (
            IOobject::groupName("kappa", phase)
        );
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
    limitGrad_(false),
    TRefType_(KAPPA_BY_DELTA),
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
    limitGrad_(dict.lookupOrDefault("limitGrad", false)),
    TRefType_
    (
        dict.found("TRefType")
      ? TRefTypeNames_.read(dict.lookup("TRefType"))
      : KAPPA_BY_DELTA
    ),
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

//     fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

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
    limitGrad_(psf.limitGrad_),
    TRefType_(psf.TRefType_),
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
    limitGrad_(psf.limitGrad_),
    TRefType_(psf.TRefType_),
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


    if (!returnReduce(nbrPatch.size(), sumOp<label>()))
    {
        refGrad() = Zero;
        valueFraction() = 0.0;
        refValue() = unmappedT_;
        mixedFvPatchScalarField::updateCoeffs();
        UPstream::msgType() = oldTag;
        return;
    }

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

    // Values for this patch
    scalarField& Tp = *this;
    scalarField TcOwn(this->patchInternalField());

    //- Values for neighbour patch
    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);
    scalarField TcNbr(nbrTp.patchInternalField());

    // Swap to obtain full local values of neighbour K*delta
    scalarField nbrKappaByDelta;
    scalarField nbrKappaTByDelta;
    {
        if (contactRes_ == 0.0)
        {
            nbrKappaByDelta = kappa(nbrTp)*nbrPatch.deltaCoeffs();
        }
        else
        {
            nbrKappaByDelta.setSize(this->size(), contactRes_);
        }

        nbrKappaTByDelta = samplePatch.faceInterpolate
        (
            nbrKappaByDelta*TcNbr
        );
        nbrKappaByDelta = samplePatch.faceInterpolate(nbrKappaByDelta);
    }

    scalarField pkappa(this->kappa(*this));
    scalarField kappaByDelta(pkappa*patch().deltaCoeffs());

    scalarField q(Tp.size(), 0.0);
    {
        const bool hOwn = hName_ != "none";
        const bool hNbr = hNbrName_ != "none";
        scalarField deltaT;
        if (hOwn || hNbr)
        {
            //- Difference in temperature
            deltaT = samplePatch.faceInterpolate(TcNbr) - TcOwn;
        }

        if (hOwn)
        {
            q +=
                patch().lookupPatchField<volScalarField, scalar>(hName_)
               *deltaT;
        }
        if (hNbr)
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
    }

    bool nbrRad = false;
    if (qrName_ != "none")
    {
        q += patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    if (qrNbrName_ != "none")
    {
        nbrRad = true;
        q +=
            samplePatch.faceInterpolate
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>
                (
                    qrNbrName_
                )
            );
    }

    scalarField& grad = refGrad();
    scalarField& vf = valueFraction();
    scalarField& rv = refValue();
    forAll(grad, i)
    {
        vf[i] = nbrKappaByDelta[i]/(nbrKappaByDelta[i] + kappaByDelta[i] + small);

        switch (TRefType_)
        {
            case KAPPA_BY_DELTA:
                if (nbrKappaByDelta[i] > small)
                {
                    rv[i] = (nbrKappaTByDelta[i] + q[i])/nbrKappaByDelta[i];
                }
                break;
            case MAX:
                rv[i] = max(TcNbr[i], TcOwn[i]);
                break;
            case NEIGHBOUR:
                rv[i] = TcNbr[i];
                break;
        }

        // Set gradient
        if (pkappa[i] > small)
        {
            grad[i] = q[i]/pkappa[i];
        }
        else
        {
            grad[i] = 0.0;
        }
    }

    // rv = max(samplePatch.faceInterpolate(TcNbr), TcOwn);
    // grad = 0.0;
    if (limitGrad_)
    {
        TcNbr = samplePatch.faceInterpolate(TcNbr);

        // Limit gradients based on neighbour cells and max/min
        // coupled region
        scalar gMinT = great;
        scalar gMaxT = -great;
        if (nbrRad)
        {
            gMinT = gMin(nbrTp.internalField());
            gMaxT = gMax(nbrTp.internalField());
        }

        const scalarField& dc = patch().deltaCoeffs();

        //- Make sure resulting temperature is within  physical bounds
        forAll(grad, i)
        {
            const scalar gf = 1.0 - vf[i];
            if (vf[i] > small && gf > small)
            {
                const scalar minT = min(min(TcOwn[i], TcNbr[i]), gMinT);
                const scalar maxT = max(max(TcOwn[i], TcNbr[i]), gMaxT);

                scalar minGradT =
                    ((minT - vf[i]*rv[i])/gf - TcOwn[i])*dc[i];
                scalar maxGradT =
                    ((maxT - vf[i]*rv[i])/gf - TcOwn[i])*dc[i];
                grad[i] = min(max(grad[i], minGradT), maxGradT);
            }
        }
    }

    if (cgpp.hasUnmappedFaces())
    {
        cgpp.setUnmappedFace(valueFraction(), 0.0);
        cgpp.setUnmappedFace(refGrad(), 0.0);
        cgpp.setUnmappedFace(refValue(), unmappedT_);
    }

    mixedFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(pkappa*patch().magSf()*snGrad());

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
    writeEntryIfDifferent<word>(os, "Tnbr", "T", TnbrName_);
    writeEntryIfDifferent<word>(os, "qrNbr", "none", qrNbrName_);
    writeEntryIfDifferent<word>(os, "qr", "none", qrName_);
    writeEntryIfDifferent<word>(os, "hNbr", "none", hNbrName_);
    writeEntryIfDifferent<word>(os, "h", "none", hName_);
    if (thicknessLayers_.size())
    {
        writeEntry(os, "thicknessLayers", thicknessLayers_);
    }
    if (kappaLayers_.size())
    {
        writeEntry(os, "kappaLayers", kappaLayers_);
    }
    writeEntryIfDifferent<scalar>
    (
        os,
        "unmappedT",
        constant::thermodynamic::Tstd,
        unmappedT_
    );
    writeEntryIfDifferent<scalar>
    (
        os,
        "limitGrad",
        false,
        limitGrad_
    );
    writeEntryIfDifferent<word>
    (
        os,
        "TRefType",
        TRefTypeNames_[KAPPA_BY_DELTA],
        TRefTypeNames_[TRefType_]
    );
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
