/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
02-06-2020  Jeff Heylmun    : Modified ODE system to solve radiation
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

#include "radiationODE.H"
#include "radiationModel.H"
#include "basicThermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationODE::radiationODE
(
    const radiationModel& rad,
    const fvMesh& mesh
)
:
    ODESystem(),
    rad_(rad),
    thermo_(mesh.lookupObject<basicThermoModel>("basicThermo")),
    dict_(rad.optionalSubDict("radiationODECoeffs")),
    solve_(dict_.lookupOrDefault("solveODE", false)),
    nEqns_(1),
    q_(1, 0.0),
    dqdt_(1, 0.0),
    deltaT_(mesh.nCells(), mesh.time().deltaTValue())
{
    if (solve_)
    {
        odeSolver_ = ODESolver::New(*this, dict_);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationODE::~radiationODE()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::radiationODE::derivatives
(
    const scalar time,
    const scalarField& q,
    scalarField& dqdt
) const
{
    scalar e = q[0]/max(thermo_.rho()[celli_], 1e-10);
    scalar T = thermo_.TRhoEi(thermo_.T()[celli_], e, celli_);
    dqdt = 0.0;
    dqdt[0] = min(rad_.Ru(celli_) - rad_.Rp(celli_)*pow4(T), -q_[0]);
}


void Foam::radiationODE::jacobian
(
    const scalar t,
    const scalarField& q,
    scalarField& dqdt,
    scalarSquareMatrix& J
) const
{
    scalar e = q[0]/max(thermo_.rho()[celli_], 1e-10);
    scalar T = thermo_.TRhoEi(thermo_.T()[celli_], e, celli_);
    dqdt = 0.0;
    dqdt[0] = min(rad_.Ru(celli_) - rad_.Rp(celli_)*pow4(T), -q_[0]);
    J = scalarSquareMatrix(1, 0.0);
}


Foam::scalar Foam::radiationODE::solve
(
    const scalar& deltaT,
    volScalarField& rhoE
)
{
    if (!odeSolver_.valid())
    {
        return min(deltaT_);
    }

    for(celli_ = 0; celli_ < rhoE.size(); celli_++)
    {
        q_ = 0.0;
        q_[0] = rhoE[celli_];

        scalar timeLeft = deltaT;
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            odeSolver_->solve(0, dt, q_, deltaT_[celli_]);
            timeLeft -= dt;
            deltaT_[celli_] = dt;
        }
        rhoE[celli_] = q_[0];
    }
    return min(deltaT_);
}


// ************************************************************************* //
