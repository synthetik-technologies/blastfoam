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

#include "neoHookeanElasticMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanElasticMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, neoHookeanElasticMisesPlastic, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar neoHookeanElasticMisesPlastic::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label neoHookeanElasticMisesPlastic::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar neoHookeanElasticMisesPlastic::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar neoHookeanElasticMisesPlastic::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::neoHookeanElasticMisesPlastic::curYieldStress
(
    const scalar curEpsilonPEq,    // Current equivalent plastic strain
    const scalar J                 // Current Jacobian
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*stressPlasticStrainSeries_->value(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::neoHookeanElasticMisesPlastic::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar,            // Scaled shear modulus
    const scalar J                 // Current Jacobian
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
           *curYieldStress
            (
                epsilonPEqOld + sqrtTwoOverThree_*DLambda,
                J
            );
}


void Foam::neoHookeanElasticMisesPlastic::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar J,                // Current Jacobian
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar, J);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar, J
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar,  J);

        if (i == MaxNewtonIter_)
        {
            WarningIn("neoHookeanElasticMisesPlastic::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Note: we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda,  J
        )/J;
}


Foam::tmp<Foam::volScalarField> Foam::neoHookeanElasticMisesPlastic::Ibar
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


Foam::tmp<Foam::surfaceScalarField> Foam::neoHookeanElasticMisesPlastic::Ibar
(
    const surfaceSymmTensorField& devBEbar
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

    tmp<surfaceScalarField> tIbar
    (
        new surfaceScalarField
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
    surfaceScalarField& Ibar = tIbar.ref();

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

    return tIbar;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::neoHookeanElasticMisesPlastic::neoHookeanElasticMisesPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    stressPlasticStrainSeries_(),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialYieldStress",
            dimPressure,
            0.0
        )
    ),
    sigmaYf_
    (
        IOobject
        (
            "sigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "initialYieldStress",
            dimPressure,
            0.0
        )
    ),
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    DSigmaYf_
    (
        IOobject
        (
            "DSigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonPf_
    (
        IOobject
        (
            "epsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    bEbarTrialf_
    (
        IOobject
        (
            "bEbarTrialf",
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
    bEbarf_
    (
        IOobject
        (
            "bEbarf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    DEpsilonPEq_
    (
        IOobject
        (
            "DEpsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DEpsilonPEqf_
    (
        IOobject
        (
            "DEpsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambdaf_
    (
        IOobject
        (
            "DLambdaf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEqf_
    (
        IOobject
        (
            "epsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
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
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    plasticNf_
    (
        IOobject
        (
            "plasticNf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    nonLinearPlasticity_(!dict.found("Hp")),
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
    bEbar_.oldTime();

    Info<< "    updateBEbarConsistent: " << updateBEbarConsistent_ << endl;

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        const dimensionedScalar E("E", dimPressure, dict);

        // Read the Poisson's ratio
        const dimensionedScalar nu("nu", dimless, dict);

        // Set the shear modulus
        mu_ = E/(2.0*(1.0 + nu));

        // Set the bulk modulus
        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_.read(dict);
        K_.read(dict);
    }
    else
    {
        FatalErrorInFunction
            << "Either E and nu or mu and K elastic parameters should be "
            << "specified"
            << abort(FatalError);
    }

    // Check if plasticity is a nonlinear function of plastic strain
    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
        stressPlasticStrainSeries_ =
            Function1<scalar>::New
            (
                "stressPlasticStrainSeries",
                dict
            );
        dimensionedScalar sigmaY
        (
            "sigmaY",
            dimPressure,
            stressPlasticStrainSeries_->value(0.0)
        );
        sigmaY_ == sigmaY;
        sigmaYf_ == sigmaY;
    }
    else
    {
        // Define linear plastic modulus
        Hp_ = dict.lookup<scalar>("Hp");
        dimensionedScalar sigmaY
        (
            "sigmaY",
            dimPressure,
            dict
        );
        sigmaY_ == sigmaY;
        sigmaYf_ == sigmaY;
    }

    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanElasticMisesPlastic::~neoHookeanElasticMisesPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    const volScalarField muBar(Ibar*mu_);

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*muBar*DLambda_/magSTrial));

    return volScalarField::New
    (
        "impK",
        //mesh(),
        //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
        //zeroGradientFvPatchScalarField::typeName
        scaleFactor*(4.0/3.0)*mu_ + K_
    );
}


Foam::tmp<Foam::scalarField>
Foam::neoHookeanElasticMisesPlastic::impK(const label patchi) const
{
    const symmTensorField& pbEbarTrial
    (
        bEbarTrial_.boundaryField()[patchi]
    );
    // Calculate deviatoric trial stress
    const symmTensorField sTrial(mu_.value()*dev(pbEbarTrial));

    const scalarField Ibar(tr(pbEbarTrial)/3.0);
    const scalarField muBar(Ibar*mu_.value());

    // Magnitude of the deviatoric trial stress
    const scalarField magSTrial(max(mag(sTrial), small));

    // Calculate scaling factor
    const scalarField scaleFactor
    (
        1.0
      - (2.0*muBar*DLambda_.boundaryField()[patchi]/magSTrial)
    );

    return scaleFactor*(4.0/3.0)*mu_.value() + K_.value();
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlastic::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        K_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlastic::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        K_ + (4.0/3.0)*mu_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::neoHookeanElasticMisesPlastic::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


void Foam::neoHookeanElasticMisesPlastic::correct(volSymmTensorField& sigma)
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
    const volTensorField relFbar(1.0/cbrt(relJ)*relF());

    // Update bE trial
    bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Calculate trial deviatoric stress
    const volSymmTensorField sTrial(mu_*dev(bEbarTrial_));

    const volScalarField Ibar(tr(bEbarTrial_)/3.0);
    const volScalarField muBar(Ibar*mu_);

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrial_.primitiveField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const volScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*J*sigmaY_);

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticN_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaY_.primitiveFieldRef();
    scalarField& DLambdaI = DLambda_.primitiveFieldRef();
    const scalarField& muBarI = muBar.primitiveField();
    const scalarField& JI = J.primitiveField();
    const scalarField& sigmaYI = sigmaY_.primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().primitiveField();

    // Calculate DLambda_ and plasticN_
    forAll(fTrialI, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);
        if (magS > SMALL)
        {
            plasticNI[cellI] = sTrialI[cellI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain where t is start of time-step
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaI[cellI],
                    curSigmaY,
                    epsilonPEqOldI[cellI],
                    magS,
                    muBarI[cellI],
                    JI[cellI],
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[cellI] = fTrialI[cellI]/(2*muBarI[cellI]);

                if (magHp > SMALL)
                {
                    DLambdaI[cellI] /= 1.0 + Hp_/(3*muBarI[cellI]);

                    // Update increment of yield stress
                    DSigmaYI[cellI] = DLambdaI[cellI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryFieldRef()[patchI];
        const scalarField& muBarP = muBar.boundaryField()[patchI];
        const scalarField& JP = J.boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = sTrialP[faceI]/magS;
            }

            // Calculate DEpsilonPEq
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
            }
            else
            {
                // yielding
                if (nonLinearPlasticity_)
                {
                    scalar curSigmaY = 0.0; // updated in loop below

                    // Calculate DEpsilonPEq and curSigmaY
                    newtonLoop
                    (
                        DLambdaP[faceI],
                        curSigmaY,
                        epsilonPEqOldP[faceI],
                        magS,
                        muBarP[faceI],
                        JP[faceI],
                        maxMagBE
                    );

                    // Update increment of yield stress
                    DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                }
                else
                {
                    // Plastic modulus is linear
                    DLambdaP[faceI] = fTrialP[faceI]/(2.0*muBarP[faceI]);

                    if (magHp > SMALL)
                    {
                        DLambdaP[faceI] /= 1.0 + Hp_/(3.0*muBarP[faceI]);

                        // Update increment of yield stress
                        DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                    }
                }
            }
        }
    }

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;
    DEpsilonP_ = Ibar*DLambda_*plasticN_;
    DEpsilonP_.relax();

    // Calculate deviatoric stress
    const volSymmTensorField s(sTrial - 2*mu_*DEpsilonP_);

    // Update bEbar
    if (updateBEbarConsistent_)
    {
        const volSymmTensorField devBEbar(s/mu_);
        bEbar_ = devBEbar + this->Ibar(devBEbar)*I;
    }
    else
    {
        bEbar_ = (s/mu_) + Ibar*I;
    }

    // Update hydrostatic stress (negative of pressure)
    updateSigmaHyd
    (
        0.5*K_*(sqr(J) - 1.0),
        (4.0/3.0)*mu_ + K_
    );

    // Update the Cauchy stress
    sigma = (1.0/J)*(sigmaHyd()*I + s);
}


void Foam::neoHookeanElasticMisesPlastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    const surfaceScalarField& Jf = mechanicalLaw::Jf();

    // Calculate the relative Jacobian
    const surfaceScalarField relJ(Jf/Jf.oldTime());

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    surfaceTensorField relFbar(1.0/cbrt(relJ)*relFf());

    // Update bE trial
    bEbarTrialf_ = transform(relFbar, bEbarf_.oldTime());

    // Calculate trial deviatoric stress
    surfaceSymmTensorField sTrial(mu_*dev(bEbarTrialf_));

    const surfaceScalarField Ibar(tr(bEbarTrialf_)/3.0);
    const surfaceScalarField muBar(Ibar*mu_);

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonPf_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrialf_.primitiveField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const surfaceScalarField fTrial
    (
        mag(sTrial) - sqrtTwoOverThree_*Jf*sigmaYf_
    );

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticNf_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaYf_.primitiveFieldRef();
    scalarField& DLambdaI = DLambdaf_.primitiveFieldRef();
    const scalarField& muBarI = muBar.primitiveField();
    const scalarField& JI = Jf.primitiveField();
    const scalarField& sigmaYI = sigmaYf_.primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().primitiveField();

    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrialI, faceI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[faceI]);
        if (magS > SMALL)
        {
            plasticNI[faceI] = sTrialI[faceI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[faceI] < SMALL)
        {
            // elastic
            DSigmaYI[faceI] = 0.0;
            DLambdaI[faceI] = 0.0;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain where t is start of time-step
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaI[faceI],
                    curSigmaY,
                    epsilonPEqOldI[faceI],
                    magS,
                    muBarI[faceI],
                    JI[faceI],
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[faceI] = curSigmaY - sigmaYI[faceI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[faceI] = fTrialI[faceI]/(2*muBarI[faceI]);

                if (magHp > SMALL)
                {
                    DLambdaI[faceI] /= 1.0 + Hp_/(3*muBarI[faceI]);

                    // Update increment of yield stress
                    DSigmaYI[faceI] = DLambdaI[faceI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Calculate on coupled for surface fields
        //if (!fTrial.boundaryField()[patchI].coupled())
        {
            // Take references to the boundary patch fields for efficiency
            const scalarField& fTrialP = fTrial.boundaryField()[patchI];
            const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
            symmTensorField& plasticNP = plasticNf_.boundaryFieldRef()[patchI];
            scalarField& DSigmaYP = DSigmaYf_.boundaryFieldRef()[patchI];
            scalarField& DLambdaP = DLambdaf_.boundaryFieldRef()[patchI];
            const scalarField& muBarP = muBar.boundaryField()[patchI];
            const scalarField& JP = Jf.boundaryField()[patchI];
            const scalarField& sigmaYP = sigmaYf_.boundaryField()[patchI];
            const scalarField& epsilonPEqOldP =
                epsilonPEq_.oldTime().boundaryField()[patchI];

            forAll(fTrialP, faceI)
            {
                // Calculate direction plasticN
                const scalar magS = mag(sTrialP[faceI]);
                if (magS > SMALL)
                {
                    plasticNP[faceI] = sTrialP[faceI]/magS;
                }

                // Calculate DEpsilonPEq
                if (fTrialP[faceI] < SMALL)
                {
                    // elasticity
                    DSigmaYP[faceI] = 0.0;
                    DLambdaP[faceI] = 0.0;
                }
                else
                {
                    // yielding
                    if (nonLinearPlasticity_)
                    {
                        scalar curSigmaY = 0.0; // updated in loop below

                        // Calculate DEpsilonPEq and curSigmaY
                        newtonLoop
                        (
                            DLambdaP[faceI],
                            curSigmaY,
                            epsilonPEqOldP[faceI],
                            magS,
                            muBarP[faceI],
                            JP[faceI],
                            maxMagBE
                        );

                        // Update increment of yield stress
                        DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                    }
                    else
                    {
                        // Plastic modulus is linear
                        DLambdaP[faceI] = fTrialP[faceI]/(2.0*muBarP[faceI]);

                        if (magHp > SMALL)
                        {
                            DLambdaP[faceI] /= 1.0 + Hp_/(3.0*muBarP[faceI]);

                            // Update increment of yield stress
                            DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                        }
                    }
                }
            }
        }
    }

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;
    DEpsilonPf_ = Ibar*DLambdaf_*plasticNf_;
    DEpsilonPf_.relax();

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(sTrial - 2*mu_*DEpsilonPf_);

    // Update bEbar
    if (updateBEbarConsistent_)
    {
        const surfaceSymmTensorField devBEbar(s/mu_);
        bEbarf_ = devBEbar + this->Ibar(devBEbar)*I;
    }
    else
    {
        bEbarf_ = (s/mu_) + Ibar*I;
    }

    // Update the Cauchy stress
    // Note: updayeSigmaHyd is not implemented for surface fields
    sigma = (1.0/Jf)*(0.5*K_*(sqr(Jf) - 1)*I + s);
}


Foam::scalar Foam::neoHookeanElasticMisesPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().time().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).foundObject<surfaceTensorField>("Ff")
    )
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonPf_.primitiveField()
                  - DEpsilonPf_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().primitiveField()));
    }
    else
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonP_.primitiveField()
                  - DEpsilonP_.prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().primitiveField()));
    }
}


void Foam::neoHookeanElasticMisesPlastic::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;
    sigmaYf_ += DSigmaYf_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonPEqf_ += DEpsilonPEqf_;
    epsilonP_ += DEpsilonP_;
    epsilonPf_ += DEpsilonPf_;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.primitiveField(), celli)
    {
        if (DEpsilonPEq_.primitiveField()[celli] > SMALL)
        {
            activeYield_.primitiveFieldRef()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.primitiveFieldRef()[celli] = 0.0;
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchi)
    {
        if (!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
                {
                    activeYield_.boundaryFieldRef()[patchi][facei] = 1.0;
                }
                else
                {
                    activeYield_.boundaryFieldRef()[patchi][facei] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
}


Foam::scalar Foam::neoHookeanElasticMisesPlastic::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformatio gradient: already done by updateF
    // if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    // {
    //     F() = fvc::average(relFf()) & F().oldTime();
    // }
    // else
    // {
    //     F() = relF() & F().oldTime();
    // }

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon(0.5*log(symm(F().T() & F())));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon))));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP_.primitiveField();
    const symmTensorField& plasticNI = plasticN_.primitiveField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().primitiveField();
    const scalarField& epsilonEqI = epsilonEq.primitiveField();

    // Calculate error field
    const symmTensorField DEpsilonPErrorI
    (
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL)
    );

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max time integration error = "
            << maxMagDEpsilonPErr
            << endl;

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningInFunction
                << "The error in the plastic strain is lover 50 times larger "
                << "than the desired value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


// ************************************************************************* //
