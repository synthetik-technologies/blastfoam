/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "linearElasticSolidModel.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"
#include "twoDPointCorrector.H"

#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solidModels
{
    defineTypeNameAndDebug(linearElastic, 0);

    addToRunTimeSelectionTable
    (
        solidModel,
        linearElastic,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * Protected Member Function * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidModels::linearElastic::linearElastic(fvMesh& mesh)
:
    solidModel(type(), mesh),
    mechanicalProperties_
    (
        IOobject
        (
            "mechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    planeStress_(mechanicalProperties_.lookup("planeStress")),
    E_
    (
        "E",
        uniformOrRead
        (
            mesh,
            mechanicalProperties_,
            dimensionedScalar("E", dimPressure, 0.0)
        )/thermoPtr_->rho()
    ),
    nu_
    (
        uniformOrRead
        (
            mesh,
            mechanicalProperties_,
            dimensionedScalar("nu", dimless, 0.0)
        )
    ),
    mu_
    (
        "mu",
        E_/(2.0*(1.0 + nu_))
    ),
    lambda_
    (
        "lambda",
        planeStress_
      ? nu_*E_/((1.0 + nu_)*(1.0 - nu_))
      : nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))
    ),

    thermalStress_(mechanicalProperties_.lookup("thermalStress")),

    alphaPtr_
    (
        thermalStress_
      ? uniformOrRead
        (
            mesh,
            mechanicalProperties_,
            dimensionedScalar("alpha", dimensionSet(0, 0, 0, -1, 0), 0.0)
        ).ptr()
      : nullptr
    ),

    compactNormalStress_(mechanicalProperties_.lookup("compactNormalStress")),
    sigmaD_
    (
        IOobject
        (
            "sigmaD",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mu_*twoSymm(fvc::grad(D_)) + lambda_*(I*tr(fvc::grad(D_)))
    ),
    divSigmaExp_
    (
        IOobject
        (
            "divSigmaExp",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::div(sigmaD_)
    ),
    sigmaEq_
    (
        IOobject
        (
            "sigmaEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(this->sigma_)))
    )
{
    if (compactNormalStress_)
    {
        divSigmaExp_ -=
            fvc::laplacian(2.0*mu_ + lambda_, D_, "laplacian(DD,D)");
    }
    else
    {
        divSigmaExp_ -=
            fvc::div((2.0*mu_ + lambda_)*fvc::grad(D_), "div(sigmaD)");
    }

    mesh.setFluxRequired(D_.name());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidModels::linearElastic::~linearElastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::solidModels::linearElastic::evolve()
{
    int iCorr = 0;
    scalar initialResidual = 0;

    tmp<volScalarField> threeKalphaTmp;
    const volScalarField& rho(thermoPtr_->rho());
    do
    {
        if (thermalStress_)
        {
            volScalarField Cp(thermoPtr_->Cp());
            volScalarField& alpha = alphaPtr_();
            if (planeStress_)
            {
                threeKalphaTmp = E_/(1.0 - nu_)*alpha;
            }
            else
            {
                threeKalphaTmp = E_/(1.0 - 2.0*nu_)*alpha;
            }

            volScalarField& T = thermoPtr_->T();
            fvScalarMatrix TEqn
            (
                fvm::ddt(rho, Cp, T)
             ==
                fvm::laplacian(thermoPtr_->kappa()/thermoPtr_->Cv(), T)
              + fvModels_.source(rho*Cp, T)
            );

            fvConstraints_.constrain(TEqn);

            TEqn.solve();

            fvConstraints_.constrain(T);
        }

        {
            fvVectorMatrix DEqn
            (
                fvm::d2dt2(rho, D_)
            ==
                fvm::laplacian(2.0*mu_ + lambda_, D_, "laplacian(DD,D)")
              + divSigmaExp_
              + fvModels_.d2dt2(D_)
            );

            if (thermalStress_)
            {
                DEqn += fvc::grad(threeKalphaTmp()*thermoPtr_->T());
            }

            fvConstraints_.constrain(DEqn);

            initialResidual = DEqn.solve().max().initialResidual();

            if (!compactNormalStress_)
            {
                divSigmaExp_ = fvc::div(DEqn.flux());
            }
        }

        {
            volTensorField gradD(fvc::grad(D_));
            sigmaD_ = mu_*twoSymm(gradD) + (lambda_*I)*tr(gradD);

            if (compactNormalStress_)
            {
                divSigmaExp_ = fvc::div
                (
                    sigmaD_ - (2.0*mu_ + lambda_)*gradD,
                    "div(sigmaD)"
                );
            }
            else
            {
                divSigmaExp_ += fvc::div(sigmaD_);
            }
        }

        Info<< "Solid Iter " << iCorr << ": residual=" << initialResidual << endl;

    } while (initialResidual > solutionTol_ && ++iCorr < nCorr_);

    thermoPtr_->correct();

    Info<< "max(T): " << max(thermoPtr_->T()).value()
        << ", min(T): " << min(thermoPtr_->T()).value() << endl;

    sigma_ = sigmaD_;

    if (thermalStress_)
    {
        sigma_ = sigma_ - I*(threeKalphaTmp*thermoPtr_->T());
    }

    sigmaEq_ = sqrt((3.0/2.0)*magSqr(dev(sigma_)));

    volPointInterpolation::New(mesh_).interpolateDisplacement
    (
        D_,
        pointD_
    );

    Info<< "Max sigmaEq = " << max(sigmaEq_).value()
        << endl;

    return iCorr == nCorr_;

}

// ************************************************************************* //
