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

#include "blastRadiationModel.H"
#include "blastAbsorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blastRadiationModel, 0);
    defineRunTimeSelectionTable(blastRadiationModel, T);
    defineRunTimeSelectionTable(blastRadiationModel, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::blastRadiationModel::initialise()
{
    const dictionary& radDict(*this);
    Info<<radDict<<endl;

    solverFreq_ = max(1, radDict.lookupOrDefault<label>("solverFreq", 1));

    absorptionEmission_.reset
    (
        radiationModels::blastAbsorptionEmissionModel::New(radDict, mesh_).ptr()
    );
    bAbsorptionEmission_.reset
    (
        &dynamicCast<const radiationModels::blastAbsorptionEmissionModel>
        (
            absorptionEmission_()
        )
    );

    scatter_.reset(radiationModels::scatterModel::New(radDict, mesh_).ptr());

    soot_.reset(radiationModels::sootModel::New(radDict, mesh_).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastRadiationModel::blastRadiationModel(const volScalarField& T)
:
    radiationModel(T),
    radODE_(*this, T.mesh())
{}


Foam::blastRadiationModel::blastRadiationModel
(
    const word& type,
    const volScalarField& T
)
:
    radiationModel(T),
    radODE_(*this, T.mesh())
{
    // Read radiationProperties
    this->readOpt() = IOobject::MUST_READ;
    static_cast<IOdictionary&>(*this) = IOdictionary(static_cast<const IOobject&>(*this));

    coeffs_ = subOrEmptyDict(type + "Coeffs");
    solverFreq_ = 1;
    initialise();
}


Foam::blastRadiationModel::blastRadiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    radiationModel(T),
    radODE_(*this, T.mesh())
{
    // Copy the constructing dictionary
    static_cast<dictionary&>(*this) = dict;

    coeffs_ = subOrEmptyDict(type + "Coeffs");
    solverFreq_ = 1;
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::blastRadiationModel::~blastRadiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blastRadiationModel::calcRhoE
(
    const dimensionedScalar& dt,
    const volScalarField& rhoE,
    const volScalarField& rho,
    const volScalarField& e,
    const volScalarField& Cv
)
{
    if (radODE_.solve())
    {
        Info<< "Solving radiation ODE" << endl;
        tmp<volScalarField> rhoENew(rho*e);
        volScalarField K(rhoE - rhoENew());

        radODE_.solve(dt.value(), rhoENew.ref());
        rhoENew.ref() += K;
        return rhoENew;
    }

    volScalarField T3(pow3(T_));

    volScalarField den(rho + dt*4.0*this->Rp()*T3/Cv);

    volScalarField eNew
    (
        (rhoE - dt*this->Rp()*T3*(T_ - 4.0*e/Cv))/den
    );
    eNew.ref() += dt*this->Ru()/den();
    return rho*eNew;
}


// ************************************************************************* //
