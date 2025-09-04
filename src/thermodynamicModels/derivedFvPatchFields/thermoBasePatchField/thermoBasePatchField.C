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

#include "thermoBasePatchField.H"
#include "fluidBlastThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoBasePatchField, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoBasePatchField::thermoBasePatchField
(
    const fvPatch& p,
    const word& defaultGroup
)
:
    patch_(p),
    determinePhase_(true),
    phaseName_(defaultGroup)
{}


Foam::thermoBasePatchField::thermoBasePatchField
(
    const thermoBasePatchField& tbpf,
    const fvPatch& p,
    const word& defaultGroup
)
:
    patch_(p),
    determinePhase_(tbpf.determinePhase_),
    phaseName_(tbpf.phaseName_)
{
    if
    (
        determinePhase_
     && p.boundaryMesh().mesh().foundObject<basicThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                defaultGroup
            )
        )
    )
    {
        phaseName_ = defaultGroup;
    }
}


Foam::thermoBasePatchField::thermoBasePatchField
(
    const fvPatch& p,
    const dictionary& dict,
    const word& defaultGroup
)
:
    patch_(p),
    determinePhase_(dict.lookupOrDefault("determinePhase", true)),
    phaseName_(dict.lookupOrDefault<word>("phase", defaultGroup))
{
    if
    (
        determinePhase_
     && p.boundaryMesh().mesh().foundObject<basicThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                defaultGroup
            )
        )
    )
    {
        phaseName_ = defaultGroup;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoBasePatchField::~thermoBasePatchField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fluidThermo& Foam::thermoBasePatchField::thermo() const
{
    return
        patch_.db().lookupObject<fluidThermo>
        (
            IOobject::groupName
            (
                physicalProperties::typeName,
                phaseName_
            )
        );
}


Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::psi() const
{
    const fluidThermo& fThermo = thermo();
    if (isA<fluidBlastThermo>(fThermo))
    {
        return
            fThermo.p().boundaryField()[patch_.index()]
           /fThermo.rho(patch_.index());
    }
    return fThermo.psi().boundaryField()[patch_.index()];
}


Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::psiIf() const
{
    const fluidThermo& fThermo = thermo();
    if (isA<fluidBlastThermo>(fThermo))
    {
        return
            fThermo.p().boundaryField()[patch_.index()].patchInternalField()
           /fThermo.rho(patch_.index());
    }
    return fThermo.psi().boundaryField()[patch_.index()].patchInternalField();
}


Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::gamma() const
{
    const fluidThermo& fThermo = thermo();
    return fThermo.gamma
    (
        fThermo.T().boundaryField()[patch_.index()],
        patch_.index()
    );
}

Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::gammaIf() const
{
    const fluidThermo& fThermo = thermo();
    return fThermo.gamma
    (
        fThermo.T().boundaryField()[patch_.index()].patchInternalField(),
        patch_.index()
    );
}


Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::gamma
(
    const scalarField& T
) const
{
    return thermo().gamma(T, patch_.index());
}


Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::speedOfSound() const
{
    const fluidThermo& fThermo = thermo();
    if (isA<fluidBlastThermo>(fThermo))
    {
        return
            dynamicCast<const fluidBlastThermo>
            (
                fThermo
            ).speedOfSound().boundaryField()[patch_.index()];
    }
    return sqrt(gamma()/psi());
}


Foam::tmp<Foam::scalarField> Foam::thermoBasePatchField::speedOfSoundIf() const
{
    const fluidThermo& fThermo = thermo();
    if (isA<fluidBlastThermo>(fThermo))
    {
        return
            dynamicCast<const fluidBlastThermo>
            (
                fThermo
            ).speedOfSound().boundaryField()[patch_.index()];
    }
    return sqrt(gammaIf()/psiIf());
}


void Foam::thermoBasePatchField::write(Ostream& os) const
{
    writeEntryIfDifferent<word>
    (
        os,
        "phase",
        word::null,
        phaseName_
    );
}


// ************************************************************************* //
