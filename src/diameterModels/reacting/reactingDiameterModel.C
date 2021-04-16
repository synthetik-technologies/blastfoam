/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "reactingDiameterModel.H"
#include "fvc.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(reactingDiameterModel, 0);
    addToRunTimeSelectionTable(diameterModel, reactingDiameterModel, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::reactingDiameterModel::reactingDiameterModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    diameterModel(mesh, dict, phaseName),
    rate_(reactionRate::New(dict)),
    pName_(dict.lookupOrDefault("pName", word("p"))),
    TName_(dict.lookupOrDefault("TName", word("T"))),
    dVdt_
    (
        IOobject
        (
            IOobject::groupName("dVdt", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", dimVolume/dimTime, 0.0)
    )
{
    this->lookupAndInitialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diameterModels::reactingDiameterModel::~reactingDiameterModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::reactingDiameterModel::solve()
{
    const volScalarField& p
    (
        this->d_.mesh().lookupObject<volScalarField>(pName_)
    );
    const volScalarField& T
    (
        this->d_.mesh().lookupObject<volScalarField>(TName_)
    );
    solve(p, T);
}


void Foam::diameterModels::reactingDiameterModel::solve
(
    const volScalarField& pi,
    const volScalarField& T
)
{
    volScalarField VOld(this->V());
    this->storeAndBlendOld(VOld, VOld_);

    volScalarField dOld(this->d_);
    this->storeAndBlendOld(dOld, dOld_);

    volScalarField dDdt(-2.0*rate_->k(pi, T));
    this->storeAndBlendDelta(dDdt, deltad_);

    const dimensionedScalar& dT(this->d_.time().deltaT());

    //- Change dRdt to dDdt
    this->d_ = dOld + dT*dDdt;
    this->d_.max(1e-10);

    dVdt_ = (this->V() - VOld)/dT;
    volScalarField dDdtLimited((this->d_ - dOld)/dT);

    //- Compute actual delta for the time step knowing the blended
    dDdt = this->calcDelta(dDdtLimited, deltad_);
    this->storeDelta(dDdt, deltad_);

    dVdt_ = this->calcDelta(dVdt_, deltaV_);
    this->storeDelta(dVdt_, deltaV_);

    const surfaceScalarField& phi
    (
        this->d_.mesh().lookupObject<surfaceScalarField>
        (
            IOobject::groupName("phi", this->d_.group())
        )
    );

    volScalarField deltaPhid(fvc::div(phi, d_) - d_*fvc::div(phi));
    this->storeAndBlendDelta(deltaPhid, deltaPhid_);
    this->d_ = dOld - dT*(deltaPhid - dDdt);
    this->d_.max(1e-10);
}


void Foam::diameterModels::reactingDiameterModel::clearODEFields()
{
    this->clearOld(dOld_);
    this->clearOld(VOld_);
    this->clearDelta(deltad_);
    this->clearDelta(deltaPhid_);
    this->clearDelta(deltaV_);
}


Foam::tmp<Foam::volScalarField>
Foam::diameterModels::reactingDiameterModel::dMdt() const
{
    const volScalarField& rho
    (
        this->d_.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("rho", this->d_.group())
        )
    );
    return dVdt_*rho;
}

// ************************************************************************* //
