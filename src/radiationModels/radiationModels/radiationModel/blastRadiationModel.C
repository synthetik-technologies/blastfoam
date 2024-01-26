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
#include "fvmSup.H"
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

const Foam::blastRadiationModel&
Foam::blastRadiationModel::readRadiationProperties(const word& type)
{
    radiationModel::readOpt() = IOobject::MUST_READ;
    radiationModel::read();

    coeffs_ = subOrEmptyDict(type + "Coeffs");
    solverFreq_ = 1;

    return *this;
}


void Foam::blastRadiationModel::initialise()
{
    solverFreq_ = max(1, lookupOrDefault<label>("solverFreq", 1));

    absorptionEmission_.reset
    (
        radiationModels::blastAbsorptionEmissionModel::New(*this, mesh_).ptr()
    );
    bAbsorptionEmission_.reset
    (
        &dynamicCast<const radiationModels::blastAbsorptionEmissionModel>
        (
            absorptionEmission_()
        )
    );

    scatter_.reset(radiationModels::scatterModel::New(*this, mesh_).ptr());

    soot_.reset(radiationModels::sootModel::New(*this, mesh_).ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastRadiationModel::blastRadiationModel(const volScalarField& T)
:
    radiationModel(T),
    radODE_(*this, T.mesh())
{}


Foam::blastRadiationModel::blastRadiationModel(const word& type, const volScalarField& T)
:
    radiationModel(T),
    radODE_(readRadiationProperties(type), T.mesh())
{
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
    radODE_(readRadiationProperties(type), T.mesh())
{
    radiationModel::readOpt() = IOobject::MUST_READ;
    radiationModel::read();

    initialise();
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::blastRadiationModel::~blastRadiationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvScalarMatrix> Foam::blastRadiationModel::Sh
(
    const basicThermo& thermo,
    const volScalarField& he
) const
{
    if (radODE_.solve())
    {
        tmp<fvScalarMatrix> theEqn
        (
            new fvScalarMatrix(he, dimEnergy/dimTime)
        );
        fvScalarMatrix& heEqn = theEqn.ref();
        scalarField& heSource = heEqn.source();

        scalarField heNew(he);
        tmp<volScalarField> trho(thermo.rho());
        const volScalarField& rho = trho();
        const scalarField& V = he.mesh().V();

        const scalar dt = he.mesh().time().deltaTValue();
        radODE_.solve(dt, rho, heNew);

        forAll(heSource, celli)
        {
            heSource[celli] +=
                rho[celli]*V[celli]
               *(heNew[celli] - he[celli])/dt;
        }

        return theEqn;
    }

    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    return
    (
        Ru()
      - fvm::Sp(4.0*Rp()*T3/Cpv, he)
      - Rp()*T3*(T_ - 4.0*he/Cpv)
    );
}


Foam::tmp<Foam::fvScalarMatrix> Foam::blastRadiationModel::ST
(
    const dimensionedScalar& rhoCp,
    volScalarField& T
) const
{
    return
    (
        Ru()/rhoCp
      - fvm::Sp(Rp()*pow3(T)/rhoCp, T)
    );
}


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
        volScalarField eNew(e);
        volScalarField K(rhoE - rho*e);

        radODE_.solve(dt.value(), rho, eNew.ref());
        return eNew*rho + K;
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
