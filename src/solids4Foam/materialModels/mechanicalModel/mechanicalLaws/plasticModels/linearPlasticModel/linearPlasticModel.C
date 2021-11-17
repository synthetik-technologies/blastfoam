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

#include "linearPlasticModel.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearPlasticModel, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::linearPlasticModel::hydrostaticStress
(
    const volScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        volScalarField& sigmaH(sigmaHyd());

        // Store previous iteration to allow relaxation, if needed
        sigmaH.storePrevIter();

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
                "void Foam::linearPlasticModel::"
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
            pressureSmoothingScaleFactor_
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
            fvm::Sp(1.0, sigmaH)
          - fvm::laplacian(rDAf, sigmaH, "laplacian(DA,sigmaHyd)")
          + fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaH)) & mesh().Sf())
         ==
            K_*trEpsilon
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the field
        sigmaH.relax();
        return sigmaH;
    }

    return K_*trEpsilon;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::linearPlasticModel::hydrostaticStress
(
    const surfaceScalarField& trEpsilon
)
{
    if (solvePressureEqn_)
    {
        FatalErrorIn
        (
            "void Foam::linearPlasticModel::calculateHydrostaticStress\n"
            "(\n"
            "    surfaceScalarField& sigmaHyd,\n"
            "    const surfaceScalarField& trEpsilon\n"
            ")"
        )   << "'solvePressureEqn' option only implemented for volField stress "
            << "calculation" << abort(FatalError);
    }
    // Directly calculate hydrostatic stress from displacement field
    return K_*trEpsilon;

}


Foam::tmp<Foam::volSymmTensorField> Foam::linearPlasticModel::epsilon() const
{
    // Lookup gradient of displacement
    const volTensorField& gradD =
        mesh().lookupObject<volTensorField>("grad(D)");
    const volSymmTensorField& sigma =
        mesh().lookupObject<volSymmTensorField>("sigma");

    tmp<volSymmTensorField> efTmp(symm(gradD));

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        volSymmTensorField& ef(efTmp.ref());

        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearPlasticModel::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        ef.replace
        (
            symmTensor::ZZ,
          - (nu_/E_)
           *(
               sigma.component(symmTensor::XX)
             + sigma.component(symmTensor::YY)
            )
          - (
                epsilonP_.component(symmTensor::XX)
              + epsilonP_.component(symmTensor::YY)
            )
        );
    }
    return efTmp;
}


Foam::tmp<Foam::surfaceSymmTensorField>
Foam::linearPlasticModel::epsilonf() const
{
    // Lookup gradient of displacement
    const surfaceTensorField& gradD =
        mesh().lookupObject<surfaceTensorField>("grad(D)f");
    const surfaceSymmTensorField& sigma =
        mesh().lookupObject<surfaceSymmTensorField>("sigma");

    tmp<surfaceSymmTensorField> efTmp(symm(gradD));

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        surfaceSymmTensorField& ef(efTmp.ref());

        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorIn
            (
                "void Foam::linearPlasticModel::"
                "correct(surfaceSymmTensorField& sigma)"
            )   << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }

        ef.replace
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
    return efTmp;
}


Foam::tmp<Foam::symmTensorField>
Foam::linearPlasticModel::epsilon(const label patchi) const
{
    // Lookup gradient of displacement
    const tensorField& gradD =
        mesh().lookupObject<volTensorField>
        (
            "grad(D)"
        ).boundaryField()[patchi];
    const symmTensorField& sigma =
        mesh().lookupObject<volSymmTensorField>
        (
            "sigma"
        ).boundaryField()[patchi];

    tmp<symmTensorField> te(symm(gradD));

    // For planeStress, correct strain in the out of plane direction
    if (planeStress())
    {
        symmTensorField& e(te.ref());

        if (mesh().solutionD()[vector::Z] > -1)
        {
            FatalErrorInFunction
                << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << abort(FatalError);
        }
        const symmTensorField& pepsilonP =
            epsilonP_.boundaryField()[patchi];

        e.replace
        (
            symmTensor::ZZ,
           -(nu_.value()/E_.value())
           *(
               sigma.component(symmTensor::XX)
             + sigma.component(symmTensor::YY)
            )
          - (
                pepsilonP.component(symmTensor::XX)
              + pepsilonP.component(symmTensor::YY)
            )
        );
    }
    return te;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::linearPlasticModel::linearPlasticModel
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    plasticModel(name, mesh, dict, nonLinGeom)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearPlasticModel::~linearPlasticModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::linearPlasticModel::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(epsilon()));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*mu_*DLambda_/magSTrial));

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
Foam::linearPlasticModel::impK(const label patchi) const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const symmTensorField e(dev(epsilon()));

    // Calculate deviatoric trial stress
    const symmTensorField sTrial
    (
        2.0*mu_.value()
       *(e - dev(epsilonP_.oldTime().boundaryField()[patchi]))
    );

    // Magnitude of the deviatoric trial stress
    const scalarField magSTrial(max(mag(sTrial), small));

    // Calculate scaling factor
    const scalarField scaleFactor
    (
        1.0
      - 2.0*mu_.value()
       *DLambda_.boundaryField()[patchi]
       /magSTrial
    );

    return scaleFactor*(4.0/3.0)*mu_.value() + K_.value();
}

void Foam::linearPlasticModel::correct(volSymmTensorField& sigma)
{
    volSymmTensorField eps(epsilon());

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(eps));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP_.oldTime())));

    // Calculate the yield function
    const volScalarField fTrial
    (
        mag(sTrial) - sqrtTwoOverThree_*sigmaY_.oldTime()
    );

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(eps.primitiveField())), SMALL);

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
            mu_.value(),
            maxMagBE
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

        const scalarField& psigmaYOld =
            sigmaY_.oldTime().boundaryField()[patchI];
        const scalarField& pepsilonPEqOld =
            epsilonPEq_.oldTime().boundaryField()[patchI];

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
    const volScalarField trEpsilon(tr(eps));

    // Update the stress
    sigma = hydrostaticStress(trEpsilon)*I + s;
}


void Foam::linearPlasticModel::correct(surfaceSymmTensorField& sigma)
{
    surfaceSymmTensorField eps(epsilonf());

    // Calculate deviatoric strain
    const surfaceSymmTensorField e(dev(eps));

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial
    (
        2.0*mu_*(e - dev(epsilonPf_.oldTime()))
    );

    // Calculate the yield function
    const surfaceScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*sigmaYf_);

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(eps.primitiveField())), SMALL);

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
            mu_.value(),
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
    const surfaceScalarField trEpsilon(tr(eps));

    // Update the stress
    sigma = hydrostaticStress(trEpsilon)*I + s;
}


// ************************************************************************* //
