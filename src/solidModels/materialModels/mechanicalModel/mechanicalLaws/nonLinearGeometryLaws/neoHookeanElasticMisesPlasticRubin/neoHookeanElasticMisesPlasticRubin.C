/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "neoHookeanElasticMisesPlasticRubin.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticMisesPlasticRubin, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElasticMisesPlasticRubin, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Store sqrt(2/3) as we use it often
    scalar neoHookeanElasticMisesPlasticRubin::sqrtTwoOverThree_ =
        ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::neoHookeanElasticMisesPlasticRubin::Ibar
(
    const volSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<volScalarField> tIbar
    (
        new volScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& Ibar = tIbar.ref();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.primitiveFieldRef();
    const symmTensorField devBEbarI(devBEbar.primitiveField());

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryFieldRef()[patchI];
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElasticMisesPlasticRubin::neoHookeanElasticMisesPlasticRubin
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("zero", dimPressure, 0.0),
    k_("zero", dimPressure, 0.0),
    kappa_
    (
        IOobject
        (
            "yieldStress",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dict.lookup("initialYieldStress"))
    ),
    K_(dict.lookup("hardeningModulus")),
    P_
    (
        IOobject
        (
            "P",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
    Je_
    (
        IOobject
        (
            "Je",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("1", dimless, 1.0)
    ),
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    plasticCompaction_(dict.lookup("plasticCompaction"))
{
    // Force storage of old time for adjustable time-step calculations
    bEbar_.oldTime();
    Je_.oldTime();

    // Read elastic parameters
    // The user can specify E and nu or mu and k
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E/(2.0*(1.0 + nu));

        // Set the bulk modulus
        if (planeStress())
        {
            k_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            k_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("k"))
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        k_ = dimensionedScalar(dict.lookup("k"));
    }
    else
    {
        FatalErrorIn
        (
            "neoHookeanElasticMisesPlasticRubin::"
            "neoHookeanElasticMisesPlasticRubin::()"
        )   << "Either E and nu or mu and k elastic parameters should be "
            << "specified" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElasticMisesPlasticRubin::~neoHookeanElasticMisesPlasticRubin()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlasticRubin::impK() const
{
    return volScalarField::New
    (
        "impK",
        mesh(),
        (4.0/3.0)*mu_ + k_
    );
}


Foam::tmp<Foam::scalarField>
Foam::neoHookeanElasticMisesPlasticRubin::impK(const label patchi) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            mesh().C().boundaryField()[patchi].size(),
            (4.0/3.0)*mu_.value() + k_.value()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlasticRubin::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        k_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlasticRubin::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        (4.0/3.0)*mu_ + k_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlasticRubin::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


void Foam::neoHookeanElasticMisesPlasticRubin::correct
(
    volSymmTensorField& sigma
)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    const volScalarField& J = mechanicalLaw::J();

    // Calculate the relative Jacobian
    const volScalarField relJ(J/J.oldTime());

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    const volTensorField relFbar(pow(relJ, -1.0/3.0)*relF());

    // Update bE trial
    bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Calculate the equiavlent trial stress
    // Note: tr(bEbar & bEbar) == magSqr(bEbar)
    const volSymmTensorField devBEbarTrial(dev(bEbarTrial_));
    const volScalarField sigmaEqTrial
    (
        mu_*sqrt((3.0/2.0)*tr(devBEbarTrial & devBEbarTrial))
    );

    // Store previous iteration for calculation of the residual
    lambda_.storePrevIter();

    // Calculate the plastic multiplier
    // Note: kappa is assumed here to be the Kirchhoff yield stress; we could
    // scale it by J to specify the Cauchy yield stress
    lambda_ =
        min
        (
            (
                (kappa_.oldTime()/sigmaEqTrial) + (K_/(3.0*mu_))
            )/(1.0 + (K_/(3.0*mu_))),
            1.0
        );

    // Deviatoric component of bEbar
    const volSymmTensorField devBEbar(lambda_*dev(bEbarTrial_));

    // Calculate the hardening
    kappa_ = kappa_.oldTime() + (1.0 - lambda_)*(K_/(3.0*mu_))*sigmaEqTrial;

    // Calcualte Je
    if (plasticCompaction_)
    {
        Je_ =
            relJ*Je_.oldTime()
            *exp
            (
                (3.0/2.0)*((1.0 - lambda_)/lambda_)
                *(
                    1.0 - (9.0/(tr(bEbar_)*tr(inv(bEbar_))))
                )
            );
    }
    else
    {
        //Je_ = 1.0*J();
        Je_ = relJ*Je_.oldTime();
    }

    // Reciprocal of J
    const volScalarField rJ(1.0/J);

    // Calculate the deviatoric Cauchy stress
    const volSymmTensorField devT(rJ*mu_*dev(bEbar_));

    // Update hydrostatic pressure
    // Note: P == -sigmaHyd
    P_ = rJ*0.5*k_*(1.0 - pow(Je_, 2.0));

    // Update the Cauchy stress
    sigma = devT - P_*I;

    // Update bEBar
    bEbar_ = devBEbar + Ibar(devBEbar)*I;
}


void Foam::neoHookeanElasticMisesPlasticRubin::correct
(
    surfaceSymmTensorField& sigma
)
{
    notImplemented(type() + " for surface fields");
}


Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin::residual()
{
    // Note: we remove mag(I) when we normalise the residual so that the
    // residual is like a strain residual
    return
        gMax
        (
            mag
            (
                lambda_.primitiveField()
              - lambda_.prevIter().primitiveField()
            )
        )/0.1;
}


void Foam::neoHookeanElasticMisesPlasticRubin::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;

    // Count cells actively yielding
    int numCellsYielding = 0;

    scalarField& activeYieldI = activeYield_.primitiveFieldRef();
    const scalarField& lambdaI = lambda_.primitiveField();

    forAll(activeYieldI, cellI)
    {
        if (lambdaI[cellI] < (1.0 - SMALL))
        {
            activeYieldI[cellI] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYieldI[cellI] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchI)
    {
        if (!activeYield_.boundaryField()[patchI].coupled())
        {
            scalarField& activeYieldP = activeYield_.boundaryFieldRef()[patchI];
            const scalarField& lambdaP = lambda_.boundaryField()[patchI];

            forAll(activeYieldP, faceI)
            {
                if (lambdaP[faceI] < (1.0 - SMALL))
                {
                    activeYieldP[faceI] = 1.0;
                }
                else
                {
                    activeYieldP[faceI] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;

    if (mesh().time().outputTime())
    {
        Info<< "Writing tauEq" << endl;

        // Calculate the deviatoric stress
        const volSymmTensorField s(mu_*dev(bEbar_));

        // Calculate the Kirchhoff stress
        const volSymmTensorField tau("tau", s - P_*I);

        // Calculate the equivalent Kirchhoff stress
        const volScalarField tauEq("tauEq", sqrt((3.0/2.0)*magSqr(dev(tau))));

        // Write the field
        tauEq.write();

        // Plastic component of total Jacobian
        Info<< "Writing Jp" << endl;

        const volScalarField Jp("Jp", J()/Je_);

        // Write the field
        Jp.write();
    }
}


Foam::scalar Foam::neoHookeanElasticMisesPlasticRubin::newDeltaT()
{
    return mesh().time().endTime().value();
}


// ************************************************************************* //
