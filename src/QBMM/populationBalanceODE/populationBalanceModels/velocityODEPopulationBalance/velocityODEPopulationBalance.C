/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018-2019 Alberto Passalacqua
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

#include "velocityODEPopulationBalance.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace PDFTransportModels
{
namespace populationBalanceModels
{
    defineTypeNameAndDebug(velocityODEPopulationBalance, 0);
    addToRunTimeSelectionTable
    (
        ODEPopulationBalanceModel,
        velocityODEPopulationBalance,
        dictionary
    );
}
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::velocityODEPopulationBalance
(
    const word& name,
    const dictionary& dict,
    const surfaceScalarField& phi
)
:
    velocityPDFODETransportModel(name, dict, phi.mesh(), "R"),
    ODEPopulationBalanceModel(name, dict, phi),
    odeType(phi.mesh(), dict),
    collision_(dict.lookup("collision")),
    collisionKernel_
    (
        Foam::populationBalanceSubModels::collisionKernel::New
        (
            dict.subDict("collisionKernel"),
            phi_.mesh(),
            quadrature_
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::~velocityODEPopulationBalance()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::updateImplicitMomentSource()
{
    if (!collision_)
    {
        return;
    }

    return collisionKernel_->updateFields();
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::implicitMomentSource
(
    const volVelocityMoment& moment
)
{
    if (!collision_)
    {
        return tmp<fvScalarMatrix>
        (
            new fvScalarMatrix
            (
                moment,
                moment.dimensions()*dimVolume/dimTime
            )
        );
    }

    return collisionKernel_->implicitCollisionSource(moment);
}


void Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::explicitMomentSource()
{
    if (!collision_ || collisionKernel_->implicit())
    {
        return;
    }

    return odeType::solve(quadrature_, 0);
}


void
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::updateCellMomentSource(const label celli)
{
    if (!collision_)
    {
        return;
    }

    return collisionKernel_->updateCells(celli);
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityODEPopulationBalance::cellMomentSource
(
    const labelList& momentOrder,
    const label celli,
    const velocityQuadratureApproximation&,
    const label
)
{
    return collisionKernel_->explicitCollisionSource(momentOrder, celli);
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityODEPopulationBalance::realizableCo() const
{
    return velocityPDFODETransportModel::realizableCo();
}


Foam::scalar Foam::PDFTransportModels::populationBalanceModels
::velocityODEPopulationBalance::CoNum() const
{
    return velocityPDFODETransportModel::CoNum();
}


bool
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::solveMomentSources() const
{
    return odeType::solveSources_;
}


bool
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::solveMomentOde() const
{
    return odeType::solveOde_;
}


void
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::update()
{
    velocityPDFODETransportModel::update();
}


void
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::solve()
{
    velocityPDFODETransportModel::solve();
}


void
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::postUpdate()
{
    collisionKernel_->preUpdate();
    velocityPDFODETransportModel::postUpdate();
}


void
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::clearODEFields()
{
    velocityPDFODETransportModel::clearODEFields();
}


bool
Foam::PDFTransportModels::populationBalanceModels::velocityODEPopulationBalance
::readIfModified()
{
    odeType::read
    (
        populationBalanceProperties_.subDict(type() + "Coeffs")
    );

    return true;
}


// ************************************************************************* //
