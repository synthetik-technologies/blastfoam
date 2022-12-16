/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "thermalRegionSolver.H"
#include "fvCFD.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(thermal, 0);
    addToRunTimeSelectionTable(regionSolver, thermal, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::thermal::thermal
(
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    regionSolver(mesh, regions),
    thermo_
    (
        solidBlastThermo::New
        (
            mesh_,
            IOdictionary
            (
                IOobject
                (
                    "thermophysicalProperties",
                    runTime_.constant(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            )
        )
    ),
    fvModels_(fvModels::New(mesh_)),
    fvConstraints_(fvConstraints::New(mesh_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::thermal::~thermal()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::thermal::initialiseMesh(const IterType)
{
    dynMesh_.update();
    if (mesh_.moving())
    {
        const_cast<surfaceScalarField&>(mesh_.phi()) == Zero;
    }
}


void Foam::regionSolvers::thermal::initialise()
{}

void Foam::regionSolvers::thermal::solve()
{
    const volScalarField TOld(thermo_->T());

    tmp<volScalarField> trho = thermo_->rho();
    const volScalarField& rho = trho();
    volScalarField& e = thermo_->he();

    label maxIter = mesh_.solutionDict().lookup<label>("maxIter");
    scalar tolerance = mesh_.solutionDict().lookup<scalar>("tolerance");

    label iter = 0;
    bool lastIter = false;
    bool converged = false;
    do
    {
        Info<<"Iteration " << iter << endl;
        if (converged)
        {
            lastIter = true;
        }
        fvScalarMatrix eEqn
        (
            fvm::ddt(rho, e)
          + thermo_->divq(e)
         ==
            fvModels_.source(rho, e)
        );

        if (!lastIter)
        {
            eEqn.relax();
        }

        fvConstraints_.constrain(eEqn);

        converged = eEqn.solve().initialResidual() < tolerance;

        if (!lastIter)
        {
            e.storePrevIter();
            e.relax();
        }

        fvConstraints_.constrain(e);

    } while (++iter < maxIter && !lastIter);

    if (converged)
    {
        Info<< "Converged in " << iter << " iterations" << endl;
    }
    else
    {
        Info<< "Did not converge in " << iter << " iterations" << endl;
    }

    thermo_->correct();

    Info<< "Min/max T:" << min(thermo_->T()).value() << ' '
        << max(thermo_->T()).value() << endl;

    // scalar error =
    //     sqrt
    //     (
    //         sum(magSqr(TOld - thermo_->T())).value()
    //        /returnReduce(thermo_->T().size(), sumOp<scalar>())
    //     );
    // Info<< TOld.name() << " error for region " << this->name() << ": "
    //     << error << endl;

}


Foam::scalar Foam::regionSolvers::thermal::CoNum() const
{
    tmp<volScalarField> magKappa;
    if (thermo_->isotropic())
    {
        magKappa = thermo_->kappa();
    }
    else
    {
        magKappa = mag(thermo_->Kappa());
    }

    tmp<volScalarField> tcp = thermo_->Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> trho = thermo_->rho();
    const volScalarField& rho = trho();

    surfaceScalarField kapparhoCpbyDelta
    (
        sqr(mesh_.surfaceInterpolation::deltaCoeffs())
       *fvc::interpolate(magKappa)
       /fvc::interpolate(cp*rho)
    );

    const scalar DiNum = max(kapparhoCpbyDelta).value()*runTime_.deltaTValue();
    const scalar meanDiNum =
        average(kapparhoCpbyDelta).value()*runTime_.deltaTValue();

    Info<< "Region: " << mesh_.name() << " Diffusion Number mean: " << meanDiNum
        << " max: " << DiNum << endl;

    return DiNum;
}


Foam::scalar Foam::regionSolvers::thermal::maxCo() const
{
    return great;
}
// ************************************************************************* //
