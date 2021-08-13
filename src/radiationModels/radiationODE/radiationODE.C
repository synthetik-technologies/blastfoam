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
#include "blastRadiationModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationODE::radiationODE
(
    const blastRadiationModel& rad,
    const fvMesh& mesh
)
:
    ODESystem(),
    rad_(rad),
    thermo_(mesh.lookupObject<blastThermo>(basicThermo::dictName)),
    solve_(rad_.lookupOrDefault("solveODE", false)),
    dict_
    (
        solve_
      ? rad.subDict("radiationODECoeffs")
      : rad_
    ),
    nEqns_(1),
    q_(1, 0.0),
    dqdt_(1, 0.0),
    deltaT_
    (
        IOobject
        (
            "radiation::deltaT",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        mesh.time().deltaT()
    )
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
    const label li,
    scalarField& dqdt
) const
{
    scalar e = q[0]/max(thermo_.rhoi(li), 1e-10);
    scalar T = thermo_.THEi(e, thermo_.T()[li], li);
    dqdt = rad_.Rui(li) - rad_.Rpi(li)*pow4(T);
}


void Foam::radiationODE::jacobian
(
    const scalar t,
    const scalarField& q,
    const label li,
    scalarField& dqdt,
    scalarSquareMatrix& J
) const
{
    scalar e = q[0]/max(thermo_.rhoi(li), 1e-10);
    scalar T = thermo_.THEi(e, thermo_.T()[li], li);
    dqdt = rad_.Rui(li) - rad_.Rpi(li)*pow4(T);
    J = scalarSquareMatrix
        (
            1,
            -4.0*rad_.Rpi(li)*pow3(T)/thermo_.Cvi(li)
        );
}


Foam::scalar Foam::radiationODE::solve
(
    const scalar& deltaT,
    volScalarField& rhoE
)
{
    if (!odeSolver_.valid())
    {
        return min(deltaT_).value();
    }

    forAll(rhoE, celli)
    {
        q_ = rhoE[celli];

        scalar timeLeft = deltaT;
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            odeSolver_->solve(0, dt, q_, celli, deltaT_[celli]);
            timeLeft -= dt;
            deltaT_[celli] = dt;
        }
        rhoE[celli] = q_[0];
    }
    return min(deltaT_).value();
}


// ************************************************************************* //
