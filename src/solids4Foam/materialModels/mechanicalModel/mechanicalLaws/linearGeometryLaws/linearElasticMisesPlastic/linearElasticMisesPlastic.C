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

#include "linearElasticMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, linearElasticMisesPlastic, linGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar linearElasticMisesPlastic::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label linearElasticMisesPlastic::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar linearElasticMisesPlastic::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar linearElasticMisesPlastic::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::linearElasticMisesPlastic::updatePlasticity
(
    symmTensor& plasticN,          // Plastic return direction
    scalar& DLambda,               // Plastic multiplier increment
    scalar& DSigmaY,               // Increment of yield stress
    scalar& sigmaY,                // Yield stress
    const scalar sigmaYOld,        // Yield stress old time
    const scalar fTrial,           // Trial yield function
    const symmTensor& sTrial,      // Trial deviatoric stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagBE          // Max strain increment magnitude
) const
{
    // Calculate DLambda/DEpsilonPEq
    if (fTrial < SMALL)
    {
        // Elasticity
        plasticN = symmTensor(I);
        DLambda = 0.0;
        DSigmaY = 0.0;
        sigmaY = sigmaYOld;
    }
    else
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrial);
        if (magS > SMALL)
        {
            plasticN = sTrial/magS;
        }
        else
        {
            // Deviatoric stress is zero so plasticN value does not matter, but
            // we will set it to the identity
            plasticN = symmTensor(I);
        }

        if (nonLinearPlasticity_)
        {
            // Update plastic multiplier (DLambda) and current yield stress
            // (sigmaY)
            newtonLoop
            (
                DLambda,
                sigmaY,
                epsilonPEqOld,
                magS,
                mu_.value(),
                maxMagBE
            );

            // Update increment of yield stress
            DSigmaY = sigmaY - sigmaYOld;
        }
        else
        {
            // Update DLambda
            DLambda = fTrial/(2*mu_.value());

            // If the isotropic linear modulus is non-zero
            if (mag(Hp_) > SMALL)
            {
                DLambda /= 1.0 + Hp_/(3*mu_.value());

                // Update increment of yield stress
                DSigmaY = DLambda*Hp_;

                // Update yield stress
                sigmaY = sigmaYOld + DSigmaY;
            }
        }
    }
}


Foam::scalar Foam::linearElasticMisesPlastic::curYieldStress
(
    const scalar curEpsilonPEq    // Current equivalent plastic strain
) const
{
    return stressPlasticStrainSeries_.value(max(curEpsilonPEq, small));
}


Foam::scalar Foam::linearElasticMisesPlastic::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar            // Scaled shear modulus
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
                epsilonPEqOld + sqrtTwoOverThree_*DLambda
            );
}


void Foam::linearElasticMisesPlastic::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar);
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
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar);

        if (i == MaxNewtonIter_)
        {
            WarningIn("linearElasticMisesPlastic::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda
        );
}


void Foam::linearElasticMisesPlastic::calculateHydrostaticStress
(
    volScalarField& sigmaHyd,
    const volScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        // Store previous iteration to allow relaxation, if needed
        sigmaHyd.storePrevIter();

        // Lookup the momentum equation inverse diagonal field
        const volScalarField* ADPtr = NULL;
        if (mesh().foundObject<volScalarField>("DEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DEqnA");
        }
        else if (mesh().foundObject<volScalarField>("DDEqnA"))
        {
            ADPtr = &mesh().lookupObject<volScalarField>("DDEqnA");
        }
        else
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "calculateHydrostaticStress\n"
                "(\n"
                "    volScalarField& sigmaHyd,\n"
                "    const volScalarField& trEpsilon\n"
                ")"
            )   << "Cannot find the DEqnA or DDEqnA field: this should be "
                << "stored in the solidModel" << abort(FatalError);
        }
        const volScalarField& AD = *ADPtr;

        // Pressure diffusivity field multiple by (4.0/3.0)*mu + K, which is
        // equivalent to (2*mu + lambda)
        // Note: we can scale this coefficient by pressureSmoothingCoeff to
        // provide greater smoothing
        const surfaceScalarField rDAf
        (
            "rDAf",
            pressureSmoothingCoeff_
           *fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        // Solve pressure laplacian
        // Note: the the laplacian and div terms combine to produce a
        // third-order smoothing/dissipation term
        // If we only used the laplacian term then the smoothing/dissipation
        // would be second-order.
        // It would be interesting to see how this compares to the JST 2nd/4th
        // order dissipation term
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(1.0, sigmaHyd)
          - fvm::laplacian(rDAf, sigmaHyd, "laplacian(DA,sigmaHyd)")
          + fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaHyd)) & mesh().Sf())
         ==
            K_*trEpsilon
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the field
        sigmaHyd.relax();
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = K_*trEpsilon;
    }
}


void Foam::linearElasticMisesPlastic::calculateHydrostaticStress
(
    surfaceScalarField& sigmaHyd,
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorIn
        (
            "void Foam::linearElasticMisesPlastic::calculateHydrostaticStress\n"
            "(\n"
            "    surfaceScalarField& sigmaHyd,\n"
            "    const surfaceScalarField& trEpsilon\n"
            ")"
        )   << "'solvePressureEqn' option only implemented for volField stress "
            << "calculation" << abort(FatalError);
    }
    else
    {
        // Directly calculate hydrostatic stress from displacement field
        sigmaHyd = K_*trEpsilon;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearElasticMisesPlastic::linearElasticMisesPlastic
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
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    stressPlasticStrainSeries_("stressPlasticStrainSeries", dict),
    solvePressureEqn_(dict.lookup("solvePressureEqn")),
    pressureSmoothingCoeff_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingCoeff", 1.0)
    ),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaHydf_
    (
        IOobject
        (
            "sigmaHydf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0)
    ),
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
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_.value(0)
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
           "initialYieldStress", dimPressure, stressPlasticStrainSeries_.value(0)
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonf_
    (
        IOobject
        (
            "epsilonf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
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
    nonLinearPlasticity_(stressPlasticStrainSeries_.x()().size() > 2),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old-time fields
    epsilon_.oldTime();
    epsilonP_.oldTime();
    epsilonPf_.oldTime();
    epsilonPEq_.oldTime();
    epsilonPEqf_.oldTime();
    plasticN_.oldTime();
    sigmaY_.oldTime();
    sigmaYf_.oldTime();

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        E_.value() = dict.lookup<scalar>("E");

        // Read the Poisson's ratio
        nu_.value() = dict.lookup<scalar>("nu");

        // Set the shear modulus
        mu_ = E_/(2.0*(1.0 + nu_));

        // Set the bulk modulus
        K_ = (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_;
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        // Read shear modulus
        mu_.value() = dict.lookup<scalar>("mu");

        // Read bulk modulus
        K_.value() = dict.lookup<scalar>("K");

        // Calculate Young's modulus
        E_ = 9.0*K_*mu_/(3.0*K_ + mu_);

        // Calculate Poisson's ratio
        nu_ = (3.0*K_ - 2.0*mu_)/(2.0*(3.0*K_ + mu_));
    }
    else
    {
        FatalErrorIn
        (
            "linearElasticMisesPlastic::linearElasticMisesPlastic::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Check if plasticity is a nonlinear function of plastic strain
    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
    }
    else
    {
        if (stressPlasticStrainSeries_.x()().size() == 1)
        {
            Info<< "    Perfect Plasticity" << endl;
        }
        else
        {
            Info<< "    Plasticity is linear" << endl;

            // Define linear plastic modulus
            Hp_ =
                (
                    stressPlasticStrainSeries_.y()()[1]
                  - stressPlasticStrainSeries_.y()()[0]
                )
               /(
                    stressPlasticStrainSeries_.x()()[1]
                  - stressPlasticStrainSeries_.x()()[0]
                );
        }
    }

    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingCoeff: " << pressureSmoothingCoeff_
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticMisesPlastic::~linearElasticMisesPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::linearElasticMisesPlastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(epsilon_));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*mu_*DLambda_/magSTrial));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::linearElasticMisesPlastic::correct(volSymmTensorField& sigma)
{
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        epsilon_ = epsilon_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        epsilon_ = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilon_.replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
          - (
                epsilonP_.component(symmTensor::XX)
              + epsilonP_.component(symmTensor::YY)
            )
        );
    }

    // Calculate deviatoric strain
    const volSymmTensorField e (dev(epsilon_));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Calculate the yield function
    const volScalarField fTrial
    (
        mag(sTrial) - sqrtTwoOverThree_*sigmaY_.oldTime()
    );

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilon_.primitiveField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticN_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaY_.primitiveFieldRef();
    scalarField& DLambdaI = DLambda_.primitiveFieldRef();
    scalarField& sigmaYI = sigmaY_.primitiveFieldRef();
    const scalarField& sigmaYOldI = sigmaY_.oldTime().primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().primitiveField();

    forAll(fTrialI, cellI)
    {
        // Update plasticN, DLambda, DSigmaY and sigmaY for this cell
        updatePlasticity
        (
            plasticNI[cellI],
            DLambdaI[cellI],
            DSigmaYI[cellI],
            sigmaYI[cellI],
            sigmaYOldI[cellI],
            fTrialI[cellI],
            sTrialI[cellI],
            epsilonPEqOldI[cellI],
            mu_.value(),
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaY_.boundaryFieldRef()[patchI];

        const scalarField& sigmaYOldP =
            sigmaY_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            updatePlasticity
            (
                plasticNP[faceI],
                DLambdaP[faceI],
                DSigmaYP[faceI],
                sigmaYP[faceI],
                sigmaYOldP[faceI],
                fTrialP[faceI],
                sTrialP[faceI],
                epsilonPEqOldP[faceI],
                mu_.value(),
                maxMagBE
            );
        }
    }

    // Update DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;

    // Store previous iteration for residual calculation
    DEpsilonP_.storePrevIter();

    // Update DEpsilonP
    DEpsilonP_ = DLambda_*plasticN_;

    // Update total plastic strain
    epsilonP_ = epsilonP_.oldTime() + DEpsilonP_;

    // Update equivalent total plastic strain
    epsilonPEq_ = epsilonPEq_.oldTime() + DEpsilonPEq_;

    // Calculate deviatoric stress
    const volSymmTensorField s(sTrial - 2*mu_*DEpsilonP_);

    // Calculate the hydrostatic pressure
    const volScalarField trEpsilon(tr(epsilon_));
    calculateHydrostaticStress(sigmaHyd_, trEpsilon);

    // Update the stress
    sigma = sigmaHyd_*I + s;
}


void Foam::linearElasticMisesPlastic::correct(surfaceSymmTensorField& sigma)
{
    // Calculate total strain
    if (incremental())
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        epsilonf_ = epsilonf_.oldTime() + symm(gradDD);
    }
    else
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        epsilonf_ = symm(gradD);
    }

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearElasticMisesPlastic::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        epsilonf_.replace
        (
            symmTensor::ZZ,
           -(nu_/E_)
           *(sigma.component(symmTensor::XX) + sigma.component(symmTensor::YY))
          - (
                epsilonPf_.component(symmTensor::XX)
              + epsilonPf_.component(symmTensor::YY)
            )
        );
    }

    // Calculate deviatoric strain
    const surfaceSymmTensorField e(dev(epsilonf_));

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial
    (
        2.0*mu_*(e - dev(epsilonPf_.oldTime()))
    );

    // Calculate the yield function
    const surfaceScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*sigmaYf_);

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(epsilonf_.primitiveField())), SMALL);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.primitiveField();
    const symmTensorField& sTrialI = sTrial.primitiveField();
    symmTensorField& plasticNI = plasticNf_.primitiveFieldRef();
    scalarField& DSigmaYI = DSigmaYf_.primitiveFieldRef();
    scalarField& DLambdaI = DLambdaf_.primitiveFieldRef();
    scalarField& sigmaYI = sigmaYf_.primitiveFieldRef();
    const scalarField& sigmaYOldI = sigmaYf_.oldTime().primitiveField();
    const scalarField& epsilonPEqOldI = epsilonPEqf_.oldTime().primitiveField();

    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrialI, faceI)
    {
        // Update plasticN, DLambda, DSigmaY and sigmaY for this face
        updatePlasticity
        (
            plasticNI[faceI],
            DLambdaI[faceI],
            DSigmaYI[faceI],
            sigmaYI[faceI],
            sigmaYOldI[faceI],
            fTrialI[faceI],
            sTrialI[faceI],
            epsilonPEqOldI[faceI],
            mu_.value(),
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticNf_.boundaryFieldRef()[patchI];
        scalarField& DSigmaYP = DSigmaYf_.boundaryFieldRef()[patchI];
        scalarField& DLambdaP = DLambdaf_.boundaryFieldRef()[patchI];
        scalarField& sigmaYP = sigmaYf_.boundaryFieldRef()[patchI];
        const scalarField& sigmaYOldP =
            sigmaYf_.oldTime().boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEqf_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            updatePlasticity
            (
                plasticNP[faceI],
                DLambdaP[faceI],
                DSigmaYP[faceI],
                sigmaYP[faceI],
                sigmaYOldP[faceI],
                fTrialP[faceI],
                sTrialP[faceI],
                epsilonPEqOldP[faceI],
                mu_.value(),
                maxMagBE
            );
        }
    }

    // Update DEpsilonPEq
    DEpsilonPEqf_ = sqrtTwoOverThree_*DLambdaf_;

    // Store previous iteration for residual calculation
    DEpsilonPf_.storePrevIter();

    // Update DEpsilonP
    DEpsilonPf_ = DLambdaf_*plasticNf_;

    // Update total plastic strain
    epsilonPf_ = epsilonPf_.oldTime() + DEpsilonPf_;

    // Update equivalent total plastic strain
    epsilonPEqf_ = epsilonPEqf_.oldTime() + DEpsilonPEqf_;

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(sTrial - 2*mu_*DEpsilonPf_);

    // Calculate the hydrostatic pressure directly from the displacement
    // field
    const surfaceScalarField trEpsilon(tr(epsilonf_));
    calculateHydrostaticStress(sigmaHydf_, trEpsilon);

    // Update the stress
    sigma = sigmaHydf_*I + s;
}


Foam::scalar Foam::linearElasticMisesPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().foundObject<surfaceTensorField>("grad(D)f")
     || mesh().foundObject<surfaceTensorField>("grad(DD)f")
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


void Foam::linearElasticMisesPlastic::updateTotalFields()
{
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

    const int nTotalCells = returnReduce(mesh().nCells(), sumOp<int>());

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << nl
        << "    " << numCellsYielding << " cells ("
        << 100.0*scalar(numCellsYielding)/scalar(nTotalCells)
        << "% of the cells in this material) are actively yielding"
        << nl << endl;

    // Write out magnitude of plastic strain
    // if (mesh().time().outputTime())
    // {
    //     volScalarField epsilonPMag
    //     (
    //         IOobject
    //         (
    //             "epsilonPMag",
    //             mesh().time().timeName(),
    //             mesh(),
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         sqrt((2.0/3.0)*magSqr(dev(epsilonP_)))
    //     );

    //     epsilonPMag.write();
    // }
}


Foam::scalar Foam::linearElasticMisesPlastic::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the Euler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(epsilon_))));

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
            WarningIn
            (
                "Foam::scalar Foam::linearElasticMisesPlastic::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is over 50 times larger "
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
