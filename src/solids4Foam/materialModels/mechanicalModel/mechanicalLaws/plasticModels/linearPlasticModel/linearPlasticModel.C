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
        volScalarField& sigmaH(sigmaHydRef());

        // Store previous iteration to allow relaxation, if needed
        sigmaHydRef().storePrevIter();

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
            FatalErrorInFunction
                << "Cannot find the DEqnA or DDEqnA field: this should be "
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
        FatalErrorInFunction
            << "'solvePressureEqn' option only implemented for volField stress "
            << "calculation" << abort(FatalError);
    }
    // Directly calculate hydrostatic stress from displacement field
    return K_*trEpsilon;

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
    const volSymmTensorField e(dev(this->epsilon()));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP().oldTime())));

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial
    (
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL))
    );

    // Calculate scaling factor
    const volScalarField scaleFactor(1.0 - (2.0*mu_*DLambda()/magSTrial));

    return volScalarField::New
    (
        "impK",
        scaleFactor*(4.0/3.0)*mu_ + K_
    );
}


Foam::tmp<Foam::scalarField>
Foam::linearPlasticModel::impK(const label patchi) const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric strain
    const symmTensorField e(dev(this->epsilon().boundaryField()[patchi]));

    // Calculate deviatoric trial stress
    const symmTensorField sTrial
    (
        2.0*mu_.value()
       *(e - dev(epsilonP().oldTime().boundaryField()[patchi]))
    );

    // Magnitude of the deviatoric trial stress
    const scalarField magSTrial(max(mag(sTrial), small));

    // Calculate scaling factor
    const scalarField scaleFactor
    (
        1.0
      - 2.0*mu_.value()
       *DLambda().boundaryField()[patchi]
       /magSTrial
    );

    return scaleFactor*(4.0/3.0)*mu_.value() + K_.value();
}

void Foam::linearPlasticModel::correct
(
    volSymmTensorField& sigma,
    const bool needUpdate
)
{
    if (needUpdate)
    {
        this->updateEpsilon(epsilonRef(), nu_/E_, sigma, epsilonP());
    }

    const volSymmTensorField eps(this->epsilon());

    // Calculate deviatoric strain
    const volSymmTensorField e(dev(eps));

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial(2.0*mu_*(e - dev(epsilonP().oldTime())));

    // Calculate the yield function
    const volScalarField fTrial
    (
        mag(sTrial) - sqrtTwoOverThree_*sigmaY().oldTime()
    );

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(eps.primitiveField())), SMALL);

    // Take references to the old fields
    volSymmTensorField& plasticN = this->plasticNRef();
    volScalarField& DLambda = this->DLambdaRef();
    volScalarField& sigmaY = this->sigmaYRef();
    const volScalarField& sigmaYOld = sigmaY.oldTime();
    const volScalarField& epsilonPEqOld = epsilonPEq().oldTime();

    forAll(fTrial, cellI)
    {
        setCellValues(cellI);

        // Update plasticN, DLambda, DSigmaY and sigmaY for this cell
        this->updatePlasticity
        (
            plasticN[cellI],
            DLambda[cellI],
            sigmaY[cellI],
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
        symmTensorField& pN = plasticN.boundaryFieldRef()[patchI];
        scalarField& pDLambda = DLambda.boundaryFieldRef()[patchI];
        scalarField& psigmaY = sigmaY.boundaryFieldRef()[patchI];

        const scalarField& psigmaYOld =
            sigmaY.oldTime().boundaryField()[patchI];
        const scalarField& pepsilonPEqOld =
            epsilonPEqOld.boundaryField()[patchI];

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
    DEpsilonPEqRef() = sqrtTwoOverThree_*DLambda;

    // Store previous iteration for residual calculation
    DEpsilonPRef().storePrevIter();

    // Update DEpsilonP
    DEpsilonPRef() = DLambda*plasticN;

    // Update total plastic strain
    epsilonPRef() = epsilonP().oldTime() + DEpsilonP();

    // Update equivalent total plastic strain
    epsilonPEqRef() = epsilonPEq().oldTime() + DEpsilonPEq();

    // Calculate deviatoric stress
    const volSymmTensorField s(sTrial - 2.0*mu_*DEpsilonP());

    // Calculate the hydrostatic pressure
    const volScalarField trEpsilon(tr(eps));

    // Update the stress
    sigma = hydrostaticStress(trEpsilon)*I + s;
    sigma.correctBoundaryConditions();
}


void Foam::linearPlasticModel::correct
(
    surfaceSymmTensorField& sigma,
    const bool needUpdate
)
{
    if (needUpdate)
    {
        this->updateEpsilon(this->epsilonfRef(), nu_/E_, sigma, epsilonPf());
    }

    const surfaceSymmTensorField& eps(this->epsilonf());

    // Calculate deviatoric strain
    const surfaceSymmTensorField e(dev(eps));

    // Calculate deviatoric trial stress
    const surfaceSymmTensorField sTrial
    (
        2.0*mu_*(e - dev(epsilonPf().oldTime()))
    );

    // Calculate the yield function
    const surfaceScalarField fTrial(mag(sTrial) - sqrtTwoOverThree_*sigmaYf());

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(eps.primitiveField())), SMALL);

    // Take references to the old fields for efficiency
    surfaceSymmTensorField& plasticNf = this->plasticNfRef();
    surfaceScalarField& DLambdaf = this->DLambdafRef();
    surfaceScalarField& sigmaYf = this->sigmaYfRef();
    const surfaceScalarField& sigmaYfOld = sigmaYf.oldTime();
    const surfaceScalarField& epsilonPEqfOld = epsilonPEqf().oldTime();

    // Calculate DLambdaf_ and plasticNf_
    forAll(fTrial, faceI)
    {
        setFaceValues(faceI);

        // Update plasticN, DLambda, DSigmaY and sigmaY for this face
        updatePlasticity
        (
            plasticNf[faceI],
            DLambdaf[faceI],
            sigmaYf[faceI],
            sigmaYfOld[faceI],
            fTrial[faceI],
            sTrial[faceI],
            epsilonPEqfOld[faceI],
            mu_.value(),
            maxMagBE
        );
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& pfTrial = fTrial.boundaryField()[patchI];
        const symmTensorField& psTrial = sTrial.boundaryField()[patchI];
        symmTensorField& pN = plasticNf.boundaryFieldRef()[patchI];
        scalarField& pDLambda = DLambdaf.boundaryFieldRef()[patchI];
        scalarField& psigmaY = sigmaYf.boundaryFieldRef()[patchI];
        const scalarField& psigmaYOld = sigmaYfOld.boundaryField()[patchI];
        const scalarField& pepsilonPEqOld =
            epsilonPEqfOld.boundaryField()[patchI];


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
    DEpsilonPEqfRef() = sqrtTwoOverThree_*DLambdaf;

    // Store previous iteration for residual calculation
    DEpsilonPfRef().storePrevIter();

    // Update DEpsilonP
    DEpsilonPfRef() = DLambdaf*plasticNf;

    // Update total plastic strain
    epsilonPfRef() = epsilonPf().oldTime() + DEpsilonPf();

    // Update equivalent total plastic strain
    epsilonPEqfRef() = epsilonPEqf().oldTime() + DEpsilonPEqf();

    // Calculate deviatoric stress
    const surfaceSymmTensorField s(sTrial - 2*mu_*DEpsilonPf());

    // Calculate the hydrostatic pressure directly from the displacement
    // field
    const surfaceScalarField trEpsilon(tr(eps));

    // Update the stress
    sigma = hydrostaticStress(trEpsilon)*I + s;
}


// ************************************************************************* //
