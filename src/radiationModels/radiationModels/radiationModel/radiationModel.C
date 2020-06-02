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

#include "radiationModel.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(radiationModel, 0);
    defineRunTimeSelectionTable(radiationModel, T);
    defineRunTimeSelectionTable(radiationModel, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::IOobject Foam::radiationModel::createIOobject(const fvMesh& mesh) const
{
    IOobject io
    (
        "radiationProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (io.typeHeaderOk<IOdictionary>(true))
    {
        io.readOpt() = IOobject::MUST_READ_IF_MODIFIED;
        return io;
    }
    else
    {
        io.readOpt() = IOobject::NO_READ;
        return io;
    }
}


void Foam::radiationModel::initialise()
{
    solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    absorptionEmission_.reset
    (
        radiationModels::absorptionEmissionModel::New(*this, mesh_).ptr()
    );

    scatter_.reset(radiationModels::scatterModel::New(*this, mesh_).ptr());

    soot_.reset(radiationModels::sootModel::New(*this, mesh_).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModel::radiationModel(const volScalarField& T)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    coeffs_(dictionary::null),
    solverFreq_(0),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr),
    radODE_(*this, T.mesh())
{}


Foam::radiationModel::radiationModel(const word& type, const volScalarField& T)
:
    IOdictionary(createIOobject(T.mesh())),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr),
    radODE_(*this, T.mesh())
{
    initialise();
}


Foam::radiationModel::radiationModel
(
    const word& type,
    const dictionary& dict,
    const volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "radiationProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(T.mesh()),
    time_(T.time()),
    T_(T),
    coeffs_(subOrEmptyDict(type + "Coeffs")),
    solverFreq_(1),
    firstIter_(true),
    absorptionEmission_(nullptr),
    scatter_(nullptr),
    soot_(nullptr),
    radODE_(*this, T.mesh())
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiationModel::~radiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiationModel::correct()
{
    if (firstIter_ || (time_.timeIndex() % solverFreq_ == 0))
    {
        calculate();
        firstIter_ = false;
    }

    if (!soot_.empty())
    {
        soot_->correct();
    }
}


Foam::tmp<Foam::volScalarField> Foam::radiationModel::calcRhoE
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


bool Foam::radiationModel::read()
{
    if (regIOobject::read())
    {
        coeffs_ = subOrEmptyDict(type() + "Coeffs");

        solverFreq_ = lookupOrDefault<label>("solverFreq", 1);
        solverFreq_ = max(1, solverFreq_);

        return true;
    }
    else
    {
        return false;
    }
}


const Foam::radiationModels::absorptionEmissionModel&
Foam::radiationModel::absorptionEmission() const
{
    if (!absorptionEmission_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation absorptionEmission model, but model is "
            << "not active" << abort(FatalError);
    }

    return absorptionEmission_();
}


const Foam::radiationModels::sootModel& Foam::radiationModel::soot() const
{
    if (!soot_.valid())
    {
        FatalErrorInFunction
            << "Requested radiation sootModel model, but model is "
            << "not active" << abort(FatalError);
    }

    return soot_();
}


// ************************************************************************* //
