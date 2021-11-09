/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "temperatureCoupledImmersedBoundaryFvPatchScalarField.H"
#include "immersedBoundaryFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "thermophysicalTransportModel.H"


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::temperatureCoupledImmersedBoundaryFvPatchScalarField::
kappaImmersed() const
{
    const fvMesh& mesh = immersedMesh_;
    const label patchi = mapper_.interfaceIndex();

    const word thermoName
    (
        IOobject::groupName
        (
            basicThermo::dictName,
            phaseName_
        )
    );

    if (mesh.foundObject<fluidThermo>(thermoName))
    {
        static word ttmName(thermophysicalTransportModel::typeName);

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
            const vectorField n(mesh_.boundary()[patchi].nf());

            return n & kappa & n;
        }
        else
        {
            return thermo.kappa(patchi);
        }
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find a fluidThermo or solidThermo instance" <<nl
            << "on the immersed mesh" << nl
            << exit(FatalError);

        return scalarField::null();
    }
}


Foam::tmp<Foam::scalarField>
Foam::temperatureCoupledImmersedBoundaryFvPatchScalarField::kappa() const
{
    const word thermoName
    (
        IOobject::groupName
        (
            basicThermo::dictName,
            phaseName_
        )
    );

    tmp<volScalarField> kappa;
    if (mesh_.foundObject<fluidThermo>(thermoName))
    {
        static word ttmName(thermophysicalTransportModel::typeName);

        if (mesh_.foundObject<thermophysicalTransportModel>(ttmName))
        {
            const thermophysicalTransportModel& ttm =
                mesh_.lookupObject<thermophysicalTransportModel>(ttmName);

            kappa = ttm.kappaEff();
        }
        else
        {
            const fluidThermo& thermo =
                mesh_.lookupObject<fluidThermo>(thermoName);

            kappa = thermo.kappa();
        }
    }
    else if (mesh_.foundObject<solidThermo>(thermoName))
    {
        const solidThermo& thermo =
            mesh_.lookupObject<solidThermo>(thermoName);

        kappa = thermo.kappa();
    }
    else
    {
        FatalErrorInFunction
            << "Cannot find a fluidThermo or solidThermo instance" << nl
            << "on the main mesh" << nl
            << exit(FatalError);

        return scalarField::null();
    }
    return this->ibm_.interpolateTo(kappa().primitiveField());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledImmersedBoundaryFvPatchScalarField::
temperatureCoupledImmersedBoundaryFvPatchScalarField
(
    volScalarField& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
:
    immersedBoundaryScalarPatchField(f, dict, ibo),
    mesh_(f.mesh()),
    immersedMesh_(ibm_.immersedMesh()),
    immersedT_(immersedMesh_.lookupObject<volScalarField>("T")),
    mapper_(*(ibm_.mapper()))
{}


Foam::temperatureCoupledImmersedBoundaryFvPatchScalarField::
temperatureCoupledImmersedBoundaryFvPatchScalarField
(
    volScalarField& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo,
    const word&
)
:
    immersedBoundaryScalarPatchField(f, dict, ibo),
    mesh_(f.mesh()),
    immersedMesh_(ibm_.immersedMesh()),
    immersedT_(immersedMesh_.lookupObject<volScalarField>("T")),
    mapper_(*(ibm_.mapper()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::temperatureCoupledImmersedBoundaryFvPatchScalarField::updateCoeffs() const
{
    scalarField Tc
    (
        ibm_.patchInternalField
        (
            this->ibm_.pMesh().lookupObject<volScalarField>
            (
                IOobject::groupName
                (
                    "T",
                    this->field_.group()
                )
            )
        )
    );
    scalarField TNbr(immersedT_.boundaryField()[mapper_.interfaceIndex()]);
    TNbr = mapper_.mapImmersedToBoundary(TNbr);


    scalarField KDeltaNbr
    (
        kappaImmersed()
       *immersedMesh_.boundary()[mapper_.interfaceIndex()].deltaCoeffs()
    );
    KDeltaNbr = mapper_.mapImmersedToBoundary(KDeltaNbr);

    scalarField KDelta(kappa()*this->ibm_.deltaCoeffs());
    scalarField valueFraction
    (
        KDeltaNbr/stabilise((KDeltaNbr + KDelta), small)
    );

    values_ = valueFraction*TNbr + (1.0 - valueFraction)*Tc;
}


void Foam::temperatureCoupledImmersedBoundaryFvPatchScalarField::setValues()
{
    scalar Tavg(immersedT_.weightedAverage(immersedMesh_.V()).value());
    this->ibm_.setInternal(field_, Tavg);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeImmersedPatchTypeField
    (
        immersedBoundaryScalarPatchField,
        temperatureCoupledImmersedBoundaryFvPatchScalarField
    );
}
// ************************************************************************* //
