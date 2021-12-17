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

#include "plasticModel.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(plasticModel, 0);


// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar plasticModel::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label plasticModel::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar plasticModel::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar plasticModel::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::plasticModel::updatePlasticity
(
    symmTensor& plasticN,           // Plastic return direction
    scalar& DLambda,                // Plastic multiplier increment
    scalar& sigmaY,                 // Yield stress
    const scalar sigmaYOld,         // Yield stress old time
    const scalar fTrial,            // Trial yield function
    const symmTensor& sTrial,       // Trial deviatoric stress
    const scalar epsilonPEqOld,     // Old equivalent plastic strain
    const scalar muBar,             // Scaled shear modulus
    const scalar maxMagBE,          // Max strain increment magnitude
    const scalar J                  // Current Jacobian
) const
{
    // Calculate DLambda/DEpsilonPEq
    if (fTrial < SMALL)
    {
        // Elasticity
        plasticN = symmTensor(I);
        DLambda = 0.0;
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

        // Update plastic multiplier (DLambda) and current yield stress
        // (sigmaY)
        newtonLoop
        (
            DLambda,
            sigmaY,
            epsilonPEqOld,
            magS,
            muBar,
            maxMagBE,
            J
        );
    }
}


Foam::scalar Foam::plasticModel::yieldFunction
(
    const scalar epsilonPEqOld,     // Old equivalent plastic strain
    const scalar magSTrial,         // Deviatoric trial stress magnitude
    const scalar DLambda,           // Plastic multiplier
    const scalar muBar,             // Scaled shear modulus
    const scalar J                  // Current Jacobian
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
            epsilonPEqOld,
            epsilonPEqOld + sqrtTwoOverThree_*DLambda,
            J
        );
}


void Foam::plasticModel::newtonLoop
(
    scalar& DLambda,                // Plastic multiplier
    scalar& curSigmaY,              // Current yield stress
    const scalar epsilonPEqOld,     // Old equivalent plastic strain
    const scalar magSTrial,         // Deviatoric trial stress magnitude
    const scalar muBar,             // Scaled shear modulus
    const scalar maxMagDEpsilon,    // Max strain increment magnitude
    const scalar J                  // Current Jacobian
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction
    (
        epsilonPEqOld,
        magSTrial,
        DLambda,
        muBar,
        J
    );
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
                epsilonPEqOld,
                magSTrial,
                DLambda + finiteDiff_,
                muBar,
                J
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        if (mag(fTrialDerivative) < small)
        {
            break;
        }

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction
        (
            epsilonPEqOld,
            magSTrial,
            DLambda,
            muBar,
            J
        );

        if (i == MaxNewtonIter_)
        {
            WarningInFunction
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Note: we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY = curYieldStress
    (
        epsilonPEqOld,
        epsilonPEqOld + sqrtTwoOverThree_*DLambda,
        J
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from dictionary
Foam::plasticModel::plasticModel
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
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
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
        FatalErrorInFunction
            << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::plasticModel::~plasticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::plasticModel::bulkModulus() const
{
    return volScalarField::New
    (
        "bulkModulus",
        mesh(),
        K_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::plasticModel::elasticModulus() const
{
    return volScalarField::New
    (
        "elasticModulus",
        mesh(),
        K_ + (4.0/3.0)*mu_
    );
}


Foam::tmp<Foam::volScalarField> Foam::plasticModel::shearModulus() const
{
    return volScalarField::New
    (
        "shearModulus",
        mesh(),
        mu_
    );
}


Foam::scalar Foam::plasticModel::residual()
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
                    DEpsilonPf().primitiveField()
                  - DEpsilonPf().prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonPf().prevIter().primitiveField()));
    }
    else
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonP().primitiveField()
                  - DEpsilonP().prevIter().primitiveField()
                )
            )/gMax(SMALL + mag(DEpsilonP().prevIter().primitiveField()));
    }
}


void Foam::plasticModel::updateTotalFields()
{
    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.primitiveField(), celli)
    {
        if (DEpsilonPEq().primitiveField()[celli] > SMALL)
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
                if (DEpsilonPEq().boundaryField()[patchi][facei] > SMALL)
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

    Info<< "    Max DEpsilonPEq is " << max(DEpsilonPEq()) << nl
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


Foam::scalar Foam::plasticModel::newDeltaT()
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
    tmp<volSymmTensorField> e(calcEpsilon(epsilonP()));
    const volScalarField epsilonEq(sqrt((2.0/3.0)*magSqr(dev(e))));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP().primitiveField();
    const symmTensorField& plasticNI = plasticN().primitiveField();
    const symmTensorField& plasticNIold = plasticN().oldTime().primitiveField();
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
                << "The error in the plastic strain is over 50 times larger "
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
