/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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

#include "coupledSolidTractionFvPatchVectorField.H"
#include "uniformDimensionedFields.H"
#include "mappedPatchBase.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "compressibleMomentumTransportModel.H"
#include "incompressibleMomentumTransportModel.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::coupledSolidTractionFvPatchVectorField::viscousStress
(
    const polyMesh& mesh,
    const fvPatch& patch
) const
{
    typedef compressibleMomentumTransportModel cmpTurbModel;
    typedef incompressibleMomentumTransportModel icoTurbModel;

    if (mesh.foundObject<volSymmTensorField>("devTau"))
    {
        return
            patch.nf()
          & patch.lookupPatchField<volSymmTensorField, symmTensor>("devTau");
    }
    else if (mesh.foundObject<cmpTurbModel>(cmpTurbModel::typeName))
    {
        const cmpTurbModel& turb
        (
            mesh.lookupObject<cmpTurbModel>(cmpTurbModel::typeName)
        );

        return
            patch.nf()
          & turb.devTau()().boundaryField()[patch.index()];
    }
    else if (mesh.foundObject<volSymmTensorField>("devSigma"))
    {
        return
            (
                patch.nf()
              & patch.lookupPatchField<volSymmTensorField, symmTensor>("devSigma")
            )*rho(mesh, patch);

    }
    else if (mesh.foundObject<icoTurbModel>(icoTurbModel::typeName))
    {
        const icoTurbModel& turb
        (
            mesh.lookupObject<icoTurbModel>(icoTurbModel::typeName)
        );

        return
            (
                turb.devSigma()().boundaryField()[patch.index()]
              & patch.nf()
            )*rho(mesh, patch);
    }
    else
    {
        // For laminar flows get the velocity
        const fvPatchVectorField& Up
        (
            patch.lookupPatchField<volVectorField, vector>("U")
        );

        return mu(mesh, patch)*Up.snGrad();
    }
}


Foam::tmp<Foam::scalarField>
Foam::coupledSolidTractionFvPatchVectorField::rho
(
    const polyMesh& mesh,
    const fvPatch& patch
) const
{
    if (mesh.foundObject<volScalarField>("rho"))
    {
        return patch.lookupPatchField<volScalarField, scalar>("rho");
    }
    if (mesh.foundObject<volScalarField>("thermo:rho"))
    {
        return patch.lookupPatchField<volScalarField, scalar>("thermo:rho");
    }
    else if (mesh.foundObject<IOdictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties =
            mesh.lookupObject<IOdictionary>("transportProperties");

        return tmp<scalarField>
        (
            new scalarField
            (
                patch.size(),
                transportProperties.lookup<scalar>("rho")
            )
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(patch.size(), 1.0));
    }
}


Foam::tmp<Foam::scalarField>
Foam::coupledSolidTractionFvPatchVectorField::mu
(
    const polyMesh& mesh,
    const fvPatch& patch
) const
{
    if (mesh.foundObject<volScalarField>("thermo:mu"))
    {
        return patch.lookupPatchField<volScalarField, scalar>("thermo:mu");
    }
    else if (mesh.foundObject<volScalarField>("mu"))
    {
        return patch.lookupPatchField<volScalarField, scalar>("mu");
    }
    else if (mesh.foundObject<IOdictionary>("transportProperties"))
    {
        const IOdictionary& transportProperties =
            mesh.lookupObject<IOdictionary>("transportProperties");

        if (transportProperties.found("nu"))
        {
             return rho(mesh, patch)*transportProperties.lookup<scalar>("nu");
        }

        return tmp<scalarField>
        (
            new scalarField
            (
                patch.size(),
                transportProperties.lookup<scalar>("mu")
            )
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(patch.size(), 0.0));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledSolidTractionFvPatchVectorField::
coupledSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(p, iF),
    mpp_(p),
    pName_("p"),
    pRef_(0.0)
{}


Foam::coupledSolidTractionFvPatchVectorField::
coupledSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    solidTractionFvPatchVectorField(p, iF),
    mpp_(p),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    pRef_(dict.lookup<scalar>("pRef"))
{}


Foam::coupledSolidTractionFvPatchVectorField::
coupledSolidTractionFvPatchVectorField
(
    const coupledSolidTractionFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    solidTractionFvPatchVectorField(tdpvf, p, iF, mapper),
    mpp_(p),
    pName_(tdpvf.pName_),
    pRef_(tdpvf.pRef_)
{}


Foam::coupledSolidTractionFvPatchVectorField::
coupledSolidTractionFvPatchVectorField
(
    const coupledSolidTractionFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidTractionFvPatchVectorField(tdpvf, iF),
    mpp_(tdpvf.mpp_),
    pName_(tdpvf.pName_),
    pRef_(tdpvf.pRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledSolidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidTractionFvPatchVectorField::autoMap(m);
}


void Foam::coupledSolidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    solidTractionFvPatchVectorField::rmap(ptf, addr);
}

void Foam::coupledSolidTractionFvPatchVectorField::updateCoeffs()
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
    const polyMesh& nbrMesh = mpp_.sampleMesh();
    const label samplePatchi = mpp_.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];

    //- Lookup viscous stress and pressure fields
    vectorField nbrViscous(viscousStress(nbrMesh, nbrPatch));

    const volScalarField& volNbrP(nbrMesh.lookupObject<volScalarField>(pName_));
    scalarField nbrP(volNbrP.boundaryField()[samplePatchi] - pRef_);
    if (volNbrP.dimensions() != dimPressure)
    {
        nbrP /= rho(nbrMesh, nbrPatch);
    }

    mpp_.distribute(nbrViscous);
    mpp_.distribute(nbrP);

    this->pressure() = nbrP;

    // Flip sign since the boundary normal is opposite and the stress is dotted
    // with the neighbor boundary then mapped
    this->traction() = -nbrViscous;

    solidTractionFvPatchVectorField::updateCoeffs();
}


void Foam::coupledSolidTractionFvPatchVectorField::write(Ostream& os) const
{
    solidTractionFvPatchVectorField::write(os);
    writeEntry(os, "pName", pName_);
    writeEntry(os, "pRef", pRef_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        coupledSolidTractionFvPatchVectorField
    );
}


// ************************************************************************* //
