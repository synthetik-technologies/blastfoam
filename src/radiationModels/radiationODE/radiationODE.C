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
        solve_
      ? dict_.lookup<scalar>("initialRadDeltaT")
      : mesh.time().deltaT()
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
    const scalar rho = thermo_.cellrho(li);
    if (rho < 1e-10)
    {
        dqdt = 0.0;
    }
    else
    {
        scalar e = q[0]/rho;
        scalar T = thermo_.cellTHE(e, thermo_.T()[li], li);
        dqdt = rad_.cellRu(li) - rad_.cellRp(li)*pow4(T);
    }
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
    const scalar rho = thermo_.cellrho(li);
    if (rho < 1e-10)
    {
        dqdt = 0.0;
        J(0, 0) = 0.0;
    }
    else
    {
        scalar e = q[0]/rho;
        scalar T = thermo_.cellTHE(e, thermo_.T()[li], li);
        dqdt = rad_.cellRu(li) - rad_.cellRp(li)*pow4(T);
        J(0, 0) = -4.0*rad_.cellRp(li)*pow3(T)/thermo_.cellCv(T, li);
    }
}


Foam::scalar Foam::radiationODE::solve
(
    const scalar& deltaT,
    const scalarField& rho,
    scalarField& e
)
{
    if (!odeSolver_.valid())
    {
        return great;
    }

    forAll(e, celli)
    {
        q_ = rho[celli]*e[celli];

        scalar timeLeft = deltaT;
        while (timeLeft > small)
        {
            scalar dt = timeLeft;
            odeSolver_->solve(0, dt, q_, celli, deltaT_[celli]);
            timeLeft -= dt;
        }
        e[celli] = q_[0]/max(rho[celli], 1e-10);
    }
    return min(deltaT_).value();
}


// ************************************************************************* //
