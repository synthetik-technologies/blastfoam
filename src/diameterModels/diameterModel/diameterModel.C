/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "diameterModel.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diameterModel, 0);
    defineRunTimeSelectionTable(diameterModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModel::diameterModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    timeIntegrationSystem
    (
        IOobject::groupName("diameterModel", phaseName),
        mesh
    ),
    d_
    (
        IOobject
        (
            IOobject::groupName("d", phaseName),
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimLength, 0.0)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModel::~diameterModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModel::requireD() const
{
    typeIOobject<volScalarField> dHeader
    (
        d_.name(),
        d_.time().name(),
        d_.mesh(),
        IOobject::MUST_READ
    );
    if (!dHeader.headerOk())
    {
        FatalErrorInFunction
            << this->type() << " diameter model requires the " << d_.name()
            << "field to be specified"
            << abort(FatalError);
    }
}

Foam::tmp<Foam::volScalarField> Foam::diameterModel::A() const
{
    return Foam::constant::mathematical::pi*sqr(d_);
}


Foam::tmp<Foam::volScalarField> Foam::diameterModel::V() const
{
    return Foam::constant::mathematical::pi*pow3(d_)/6.0;
}


void Foam::diameterModel::update()
{}


void Foam::diameterModel::solve()
{}


void Foam::diameterModel::solve(const volScalarField& p, const volScalarField& T)
{
    solve();
}


void Foam::diameterModel::postUpdate()
{}


Foam::tmp<Foam::volScalarField> Foam::diameterModel::dMdt() const
{
    return volScalarField::New
    (
        IOobject::groupName("dVdt", d_.group()),
        d_.mesh(),
        dimensionedScalar(dimMass/dimTime, 0.0)
    );
}

Foam::tmp<Foam::volScalarField> Foam::diameterModel::LDratio() const
{
    return volScalarField::New
    (
        IOobject::groupName("LDratio", d_.group()),
        d_.mesh(),
        dimensionedScalar(dimless, 1.0)
    );
}

Foam::scalar
Foam::diameterModel::LDratio(const label celli) const
{
    return 1.0;
}

// ************************************************************************* //
