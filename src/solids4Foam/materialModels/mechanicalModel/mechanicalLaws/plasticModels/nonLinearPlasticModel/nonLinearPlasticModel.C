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

#include "nonLinearPlasticModel.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::nonLinearPlastic::makeJ()
{
    if (JPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    JPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                "lawJ",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // Store the old-time
    JPtr_().oldTime();
}


Foam::volScalarField& Foam::nonLinearPlastic::J()
{
    if (JPtr_.empty())
    {
        makeJ();
    }

    return JPtr_();
}


void Foam::nonLinearPlastic::makeJf()
{
    if (JfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set"
            << abort(FatalError);
    }

    JfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                "lawJf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        )
    );

    // Store the old-time
    JfPtr_().oldTime();
}


Foam::surfaceScalarField& Foam::nonLinearPlastic::Jf()
{
    if (JfPtr_.empty())
    {
        makeJf();
    }

    return JfPtr_();
}


Foam::tmp<Foam::volScalarField> Foam::nonLinearPlastic::Ibar
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


Foam::tmp<Foam::surfaceScalarField> Foam::nonLinearPlastic::Ibar
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
                alpha1 =
                    3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 =
                    3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
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


Foam::tmp<Foam::volSymmTensorField> Foam::nonLinearPlastic::epsilon() const
{
    return 0.5*log(symm(F().T() & F()));
}


Foam::tmp<Foam::surfaceSymmTensorField>
Foam::nonLinearPlastic::epsilonf() const
{
    NotImplemented;
    return bEbarf_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::nonLinearPlastic::nonLinearPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    plasticModel(name, mesh, dict, nonLinGeom),
    JPtr_(),
    JfPtr_(),
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
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    )
{
    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonLinearPlastic::~nonLinearPlastic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::nonLinearPlastic::impK() const
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


void Foam::nonLinearPlastic::correct(volSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateF(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    J() = det(F());

    // Calculate the relative Jacobian
    const volScalarField relJ(J()/J().oldTime());

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
    const scalar maxMagBE =
        max(gMax(mag(bEbarTrial_.primitiveField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const volScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*J()*sigmaY_);

    // Take references to the old fields
    const volScalarField& sigmaYOld = sigmaY_.oldTime();
    const volScalarField& epsilonPEqOld = epsilonPEq_.oldTime();

    forAll(fTrial, cellI)
    {
        setCellValues(cellI);

        // Update plasticN, DLambda, DSigmaY and sigmaY for this cell
        this->updatePlasticity
        (
            plasticN_[cellI],
            DLambda_[cellI],
            sigmaY_[cellI],
            sigmaYOld[cellI],
            fTrial[cellI],
            sTrial[cellI],
            epsilonPEqOld[cellI],
            muBar[cellI],
            maxMagBE,
            J()[cellI]
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& pfTrial = fTrial.boundaryField()[patchI];
        const symmTensorField& psTrial = sTrial.boundaryField()[patchI];
        symmTensorField& pN = plasticN_.boundaryFieldRef()[patchI];
        scalarField& pDLambda = DLambda_.boundaryFieldRef()[patchI];
        scalarField& psigmaY = sigmaY_.boundaryFieldRef()[patchI];
        const scalarField& pmuBar = muBar.boundaryField()[patchI];

        const scalarField& psigmaYOld =
            sigmaY_.oldTime().boundaryField()[patchI];
        const scalarField& pepsilonPEqOld =
            epsilonPEq_.oldTime().boundaryField()[patchI];
        const scalarField& pJ = J().boundaryField()[patchI];

        forAll(pfTrial, faceI)
        {
            setVolPatchFaceValues(patchI, faceI);

            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            this->updatePlasticity
            (
                pN[faceI],
                pDLambda[faceI],
                psigmaY[faceI],
                psigmaYOld[faceI],
                pfTrial[faceI],
                psTrial[faceI],
                pepsilonPEqOld[faceI],
                pmuBar[faceI],
                maxMagBE,
                pJ[faceI]
            );
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
        0.5*K_*(sqr(J()) - 1.0),
        (4.0/3.0)*mu_ + K_
    );

    // Update the Cauchy stress
    sigma = (1.0/J())*(sigmaHyd()*I + s);

    epsilonPEq_ += DEpsilonPEq_;
    epsilonP_ += DEpsilonP_;
}


void Foam::nonLinearPlastic::correct(surfaceSymmTensorField& sigma)
{
    // Update the deformation gradient field
    // Note: if true is returned, it means that linearised elasticity was
    // enforced by the solver via the enforceLinear switch
    if (updateFf(sigma, mu_, K_))
    {
        return;
    }

    // Update the Jacobian of the total deformation gradient
    Jf() = det(Ff());

    // Calculate the relative Jacobian
    const surfaceScalarField relJ(Jf()/Jf().oldTime());

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
        mag(sTrial) - sqrtTwoOverThree_*Jf()*sigmaYf_
    );

    // Take references to the old fields for efficiency
    const surfaceScalarField& sigmaYOld = sigmaYf_.oldTime();
    const surfaceScalarField& epsilonPEqOld = epsilonPEqf_.oldTime();

    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrial, faceI)
    {
        setFaceValues(faceI);

        // Update plasticN, DLambda, DSigmaY and sigmaY for this face
        updatePlasticity
        (
            plasticNf_[faceI],
            DLambdaf_[faceI],
            sigmaYf_[faceI],
            sigmaYOld[faceI],
            fTrial[faceI],
            sTrial[faceI],
            epsilonPEqOld[faceI],
            muBar[faceI],
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& pfTrial = fTrial.boundaryField()[patchI];
        const symmTensorField& psTrial = sTrial.boundaryField()[patchI];
        symmTensorField& pN = plasticNf_.boundaryFieldRef()[patchI];
        scalarField& pDLambda = DLambdaf_.boundaryFieldRef()[patchI];
        scalarField& psigmaY = sigmaYf_.boundaryFieldRef()[patchI];
        const scalarField& pmuBar = muBar.boundaryField()[patchI];
        const scalarField& psigmaYOld = sigmaYOld.boundaryField()[patchI];
        const scalarField& pepsilonPEqOld =
            epsilonPEqOld.boundaryField()[patchI];

        forAll(pfTrial, faceI)
        {
            setSurfacePatchFaceValues(patchI, faceI);

            // Update plasticN, DLambda, DSigmaY and sigmaY for this face
            this->updatePlasticity
            (
                pN[faceI],
                pDLambda[faceI],
                psigmaY[faceI],
                psigmaYOld[faceI],
                pfTrial[faceI],
                psTrial[faceI],
                pepsilonPEqOld[faceI],
                pmuBar[faceI],
                maxMagBE
            );
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
    sigma = (1.0/Jf())*(0.5*K_*(sqr(Jf()) - 1)*I + s);

    epsilonPEqf_ += DEpsilonPEqf_;
    epsilonPf_ += DEpsilonPf_;
}


// ************************************************************************* //
