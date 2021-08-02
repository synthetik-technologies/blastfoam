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

#include "blastTemperatureCoupledBase.H"
#include "volFields.H"
#include "fluidBlastThermo.H"
#include "solidBlastThermo.H"
#include "thermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blast::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName
)
:
    patch_(patch),
    alphaAniName_(alphaAniName)
{}


Foam::blast::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAni",word::null))
{}


Foam::blast::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const temperatureCoupledBase& base
)
:
    patch_(patch),
    alphaAniName_(base.alphaAniName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::blast::temperatureCoupledBase::kappa
(
    const fvPatchScalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    const word& phase(Tp.internalField().group());

    const word thermoName
    (
        IOobject::groupName(basicBlastThermo::typeName, phase)
    );

    if (mesh.foundObject<fluidBlastThermo>(thermoName))
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
            const fluidBlastThermo& thermo =
                mesh.lookupObject<fluidBlastThermo>(thermoName);

            return thermo.kappa(patchi);
        }
    }
    else if (mesh.foundObject<solidBlastThermo>(thermoName))
    {
        const solidBlastThermo& thermo =
            mesh.lookupObject<solidBlastThermo>(thermoName);

        if (alphaAniName_ != word::null)
        {
            const symmTensorField& alphaAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    alphaAniName_
                );

            const scalarField& rhop = thermo.rho().boundaryField()[patchi];
            const scalarField& ep = thermo.e().boundaryField()[patchi];

            const symmTensorField kappa
            (
                alphaAni*thermo.Cp(rhop, ep, Tp, patchi)
            );
            const vectorField n(patch_.nf());

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
            << "Cannot find a fluidThermo or solidThermo instance"
            << exit(FatalError);

        return scalarField::null();
    }
}


void Foam::blast::temperatureCoupledBase::write(Ostream& os) const
{
    if (alphaAniName_ != word::null)
    {
        writeEntry(os, "alphaAni", alphaAniName_);
    }
}


// ************************************************************************* //
