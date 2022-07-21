/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
31-03-2022 Synthetik Applied Technologies:  Generalized functions and added
                                            Macros for easier construction
-------------------------------------------------------------------------------
License
    This file is a derivative work of foam-extend.

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

InClass
    Foam::mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "lookupSolidModel.H"
#include "logVolFields.H"
#include "solidModel.H"
#include "fvm.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"
#include "solidSubMeshes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicalLaw, 0);
    defineRunTimeSelectionTable(mechanicalLaw, linGeomMechLaw);
    defineRunTimeSelectionTable(mechanicalLaw, nonLinGeomMechLaw);
}

// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::mechanicalLaw::makeF() const
{
    if (FPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ = makeTypeField<tensor, fvPatchField, volMesh>
    (
        "F", dimensionedTensor("I", dimless, I)
    );
}


void Foam::mechanicalLaw::makeFf() const
{
    if (FfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_ = makeTypeField<tensor, fvsPatchField, surfaceMesh>
    (
        "Ff", dimensionedTensor("I", dimless, I)
    );
}

void Foam::mechanicalLaw::makeRelF() const
{
    if (relFPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    relFPtr_ = makeTypeField<tensor, fvPatchField, volMesh>
    (
        "relF", dimensionedTensor("I", dimless, I)
    );
}


void Foam::mechanicalLaw::makeRelFf() const
{
    if (relFfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    relFfPtr_ = makeTypeField<tensor, fvsPatchField, surfaceMesh>
    (
        "relFf", dimensionedTensor("I", dimless, I)
    );
}

void Foam::mechanicalLaw::makeJ() const
{
    if (JPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    JPtr_ = makeTypeField<scalar, fvPatchField, volMesh>
    (
        "J", dimensionedScalar("J", dimless, 1.0)
    );
}


void Foam::mechanicalLaw::makeJf() const
{
    if (JfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    JfPtr_ = makeTypeField<scalar, fvsPatchField, surfaceMesh>
    (
        "Jf", dimensionedScalar("J", dimless, 1.0)
    );
}


void Foam::mechanicalLaw::makeRelJ() const
{
    if (relJPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    relJPtr_ = makeTypeField<scalar, fvPatchField, volMesh>
    (
        "relJ", dimensionedScalar("relJ", dimless, Zero)
    );
}


void Foam::mechanicalLaw::makeRelJf() const
{
    if (relJfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    relJfPtr_ = makeTypeField<scalar, fvsPatchField, surfaceMesh>
    (
        "relJf", dimensionedScalar("relJ", dimless, Zero)
    );
}

void Foam::mechanicalLaw::makeSigmaHyd() const
{
    if (sigmaHydPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    sigmaHydPtr_.set
    (
        new volScalarField
        (
            IOobject
            (
                IOobject::groupName("sigmaHyd", name_),
                mesh_.time().timeName(mesh_.time().startTime().value()),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPressure, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


void Foam::mechanicalLaw::makeSigmaHydf() const
{
    if (sigmaHydfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    sigmaHydfPtr_.set
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("sigmaHydf", name_),
                mesh_.time().timeName(mesh_.time().startTime().value()),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimPressure, 0.0)
        )
    );
}

void Foam::mechanicalLaw::makeGradSigmaHyd() const
{
    if (gradSigmaHydPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
    }

    gradSigmaHydPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                IOobject::groupName("grad(sigmaHyd)", name_),
                mesh_.time().timeName(mesh_.time().startTime().value()),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector("zero", dimPressure/dimLength, Zero),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}

// * * * * * * * * * * * * * * Protected Members * * * * * * * * * * * * * * //

bool Foam::mechanicalLaw::planeStress() const
{
    if (mesh_.foundObject<IOdictionary>("mechanicalProperties"))
    {
        return
            Switch
            (
                mesh_.lookupObject<IOdictionary>
                (
                    "mechanicalProperties"
                ).lookup("planeStress")
            );
    }
    else
    {
        // It is not straight-forward to lookup the mechanicalProperties from
        // here as we only have access to a subMesh fvMesh objectRegistry
        // We will read it here again; this switch only gets called at the start
        // of a simulation so it is not a problem
        IOdictionary mechProp
        (
            IOobject
            (
                "mechanicalProperties",
                "constant",
                mesh_.time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        return Switch(mechProp.lookup("planeStress"));
    }
}


const Foam::volTensorField& Foam::mechanicalLaw::F() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<volTensorField>("F");
    }
    if (FPtr_.empty())
    {
        makeF();
    }

    return FPtr_();
}


Foam::volTensorField& Foam::mechanicalLaw::FRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<volTensorField>("F");
    }
    if (FPtr_.empty())
    {
        makeF();
    }

    return FPtr_();
}


const Foam::surfaceTensorField& Foam::mechanicalLaw::Ff() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<surfaceTensorField>("Ff");
    }

    if (FfPtr_.empty())
    {
        makeFf();
    }

    return FfPtr_();
}


Foam::surfaceTensorField& Foam::mechanicalLaw::FfRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<surfaceTensorField>("Ff");
    }

    if (FfPtr_.empty())
    {
        makeFf();
    }

    return FfPtr_();
}


const Foam::volTensorField& Foam::mechanicalLaw::relF() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<volTensorField>("relF");
    }
    if (relFPtr_.empty())
    {
        makeRelF();
    }

    return relFPtr_();
}


Foam::volTensorField& Foam::mechanicalLaw::relFRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<volTensorField>("relF");
    }
    if (relFPtr_.empty())
    {
        makeRelF();
    }

    return relFPtr_();
}


const Foam::surfaceTensorField& Foam::mechanicalLaw::relFf() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<surfaceTensorField>("relFf");
    }
    if (relFfPtr_.empty())
    {
        makeRelFf();
    }

    return relFfPtr_();
}


Foam::surfaceTensorField& Foam::mechanicalLaw::relFfRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<surfaceTensorField>("relFf");
    }
    if (relFfPtr_.empty())
    {
        makeRelFf();
    }

    return relFfPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::J() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<volScalarField>("J");
    }

    if (JPtr_.empty())
    {
        makeJ();
    }

    return JPtr_();
}


Foam::volScalarField& Foam::mechanicalLaw::JRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<volScalarField>("J");
    }

    if (JPtr_.empty())
    {
        makeJ();
    }

    return JPtr_();
}


const Foam::surfaceScalarField&
Foam::mechanicalLaw::Jf() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<surfaceScalarField>("Jf");
    }

    if (JfPtr_.empty())
    {
        makeJf();
    }

    return JfPtr_();
}


Foam::surfaceScalarField& Foam::mechanicalLaw::JfRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<surfaceScalarField>("Jf");
    }

    if (JfPtr_.empty())
    {
        makeJf();
    }

    return JfPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::relJ() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<volScalarField>("relJ");
    }
    if (relJPtr_.empty())
    {
        makeRelJ();
    }

    return relJPtr_();
}


Foam::volScalarField& Foam::mechanicalLaw::relJRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<volScalarField>("relJ");
    }
    if (relJPtr_.empty())
    {
        makeRelF();
    }

    return relJPtr_();
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::relJf() const
{
    if (isBaseRegion())
    {
        return mesh_.lookupObject<surfaceScalarField>("relJf");
    }
    if (relJfPtr_.empty())
    {
        makeRelFf();
    }

    return relJfPtr_();
}


Foam::surfaceScalarField& Foam::mechanicalLaw::relJfRef()
{
    if (isBaseRegion())
    {
        return mesh_.lookupObjectRef<surfaceScalarField>("relJf");
    }
    if (relJfPtr_.empty())
    {
        makeRelJf();
    }

    return relJfPtr_();
}


const Foam::volScalarField& Foam::mechanicalLaw::sigmaHyd() const
{
    if (sigmaHydPtr_.empty())
    {
        makeSigmaHyd();
    }

    return sigmaHydPtr_();
}


Foam::volScalarField& Foam::mechanicalLaw::sigmaHydRef()
{
    if (sigmaHydPtr_.empty())
    {
        makeSigmaHyd();
    }

    return sigmaHydPtr_();
}


const Foam::surfaceScalarField& Foam::mechanicalLaw::sigmaHydf() const
{
    if (sigmaHydfPtr_.empty())
    {
        makeSigmaHydf();
    }

    return sigmaHydfPtr_();
}


Foam::surfaceScalarField& Foam::mechanicalLaw::sigmaHydfRef()
{
    if (sigmaHydfPtr_.empty())
    {
        makeSigmaHydf();
    }

    return sigmaHydfPtr_();
}


const Foam::volVectorField& Foam::mechanicalLaw::gradSigmaHyd() const
{
    if (gradSigmaHydPtr_.empty())
    {
        makeGradSigmaHyd();
    }

    return gradSigmaHydPtr_();
}


Foam::volVectorField& Foam::mechanicalLaw::gradSigmaHydRef()
{
    if (gradSigmaHydPtr_.empty())
    {
        makeGradSigmaHyd();
    }

    return gradSigmaHydPtr_();
}


bool Foam::mechanicalLaw::updateF
(
    volSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
{
    if (!isBaseRegion())
    {
        const fvMesh& baseMesh = this->baseMesh();
        const fvMeshSubset& subsetter =
            baseMesh.lookupObject<solidSubMeshes>
            (
                solidSubMeshes::typeName
            )[mesh().name()];
        FRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<volTensorField>("F")
        );
        relFRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<volTensorField>("relF")
        );
        relJRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<volScalarField>("relJ")
        );
        JRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<volScalarField>("J")
        );
    }

    if (enforceLinear())
    {
        WarningInFunction
            << "Material linearity enforced for stability!"
            << endl;

        // Check if the mathematical model is in total or updated
        // Lagrangian form
        if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu*symm(gradDD) + (K - (2.0/3.0)*mu)*tr(gradDD)*I;
        }
        else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
             // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
                + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;
        }
        else if (nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
        {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown nonLinGeom type: " << nonLinGeom()
                << abort(FatalError);
        }
        return true;
    }
    return false;
}


bool Foam::mechanicalLaw::updateFf
(
    surfaceSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
{
    if (!isBaseRegion())
    {
        const fvMesh& baseMesh = this->baseMesh();
        const fvMeshSubset& subsetter =
            baseMesh.lookupObject<solidSubMeshes>
            (
                solidSubMeshes::typeName
            )[mesh().name()];

        FfRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<surfaceTensorField>("Ff")
        );
        relFfRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<surfaceTensorField>("relFf")
        );
        JfRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<surfaceScalarField>("Jf")
        );
        relJfRef() = subsetter.interpolate
        (
            baseMesh.lookupObject<surfaceScalarField>("relJf")
        );
    }

    if (enforceLinear())
    {
        WarningInFunction
            << "Material linearity enforced for stability!"
            << endl;

        // Check if the mathematical model is in total or updated
        // Lagrangian form
        if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            // Lookup gradient of displacement increment
            const surfaceTensorField& gradDD =
                mesh().lookupObject<surfaceTensorField>("grad(DD)f");

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu*symm(gradDD) + (K - (2.0/3.0)*mu)*tr(gradDD)*I;
        }
        else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
             // Lookup gradient of displacement increment
            const surfaceTensorField& gradDD =
                mesh().lookupObject<surfaceTensorField>("grad(DD)f");

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
                + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;
        }
        else if (nonLinGeom() == nonLinearGeometry::LINEAR_GEOMETRY)
        {
            // Lookup gradient of displacement
            const surfaceTensorField& gradD =
                mesh().lookupObject<surfaceTensorField>("grad(D)f");

            sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown nonLinGeom type: " << nonLinGeom()
                << abort(FatalError);
        }
        return true;
    }
    return false;
}


void Foam::mechanicalLaw::updateSigmaHyd
(
    const volScalarField& sigmaHydExplicit,
    const dimensionedScalar& impK
)
{
    if (solvePressureEqn_)
    {
        SolverPerformance<scalar>::debug = 0;

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
                << "stored in the solidModel"
                << abort(FatalError);
        }
        const volScalarField& AD = *ADPtr;

        // Pressure diffusivity field
        const surfaceScalarField rDAf
        (
            "rDAf",
            pressureSmoothingScaleFactor_*fvc::interpolate
            (
                impK/AD, "interpolate(" + gradSigmaHyd().name() + ")"
            )
        );
        const dimensionedScalar one("one", dimless, 1.0);

        // Solve pressure laplacian
        // Note: the fvm and fvc laplacian terms cancel at convergence and the
        // laplacian - div(grad) term produce a smoothing/diffusion to quell
        // oscillations
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(one, sigmaHyd())
          - fvm::laplacian(rDAf, sigmaHyd(), "laplacian(rDA,sigmaHyd)")
         ==
            sigmaHydExplicit
          - fvc::div(rDAf*fvc::interpolate(gradSigmaHyd()) & mesh().Sf())
        );

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Relax the pressure field
        sigmaHydRef().relax();
    }
    else
    {
        // Explicitly calculate hydrostatic stress
        // We use 1.0* to overwritting the field IOobject attributes e.g. its
        // name and writeOpt
        sigmaHydRef() = 1.0*sigmaHydExplicit;
    }

    // Update the gradient
    gradSigmaHydRef() = fvc::grad(sigmaHyd());
}


const Foam::Switch& Foam::mechanicalLaw::enforceLinear() const
{
    // Lookup the solideModel
    const solidModel& solMod = lookupSolidModel(mesh(), baseMeshRegionName_);

    return solMod.enforceLinear();
}


bool Foam::mechanicalLaw::incremental() const
{
    // Lookup the solideModel
    const solidModel& solMod = lookupSolidModel(mesh(), baseMeshRegionName_);

    return solMod.incremental();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicalLaw::mechanicalLaw
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    name_(name),
    mesh_(mesh),
    dict_(dict),
    baseMeshRegionName_(mesh.name()),
    nonLinGeom_(nonLinGeom),
    FPtr_(),
    FfPtr_(),
    relFPtr_(),
    relFfPtr_(),
    sigmaHydPtr_(),
    gradSigmaHydPtr_(),
    useSolidDeformation_(false),
    usePlaneStress_(planeStress()),
    planeStressDir_(-1),
    nonPlaneStressDirs_(-1, -1),
    solvePressureEqn_
    (
        dict.lookupOrDefault<Switch>("solvePressureEqn", false)
    ),
    pressureSmoothingScaleFactor_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingScaleFactor", 1.0)
    )
{
    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingScaleFactor: "
            << pressureSmoothingScaleFactor_
            << endl;
    }

    if (usePlaneStress_)
    {
        const Vector<label>& solutionD(mesh_.solutionD());
        label nD = 0;
        if (solutionD[0] < 0)
        {
            nD++;
            planeStressDir_ = symmTensor::XX;
            nonPlaneStressDirs_[0] = direction(symmTensor::YY);
            nonPlaneStressDirs_[1] = direction(symmTensor::ZZ);
        }
        if (solutionD[1] < 0)
        {
            nD++;
            nonPlaneStressDirs_[0] = direction(symmTensor::XX);
            planeStressDir_ = symmTensor::YY;
            nonPlaneStressDirs_[1] = direction(symmTensor::ZZ);
        }
        if (solutionD[2] < 0)
        {
            nD++;
            nonPlaneStressDirs_[0] = direction(symmTensor::XX);
            nonPlaneStressDirs_[1] = direction(symmTensor::YY);
            planeStressDir_ = symmTensor::ZZ;
        }
        if (nD != 1)
        {
            FatalErrorInFunction
                << "For planeStress, this material law assumes one empty "
                << "direction, but " << nD << "  empty directions were found."
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * //


Foam::tmp<Foam::surfaceScalarField> Foam::mechanicalLaw::impKf() const
{
    return fvc::interpolate(impK());
}


void Foam::mechanicalLaw::correct(surfaceSymmTensorField&)
{
    notImplemented
    (
        type() + "::correct(surfaceSymmTensorField&)\n"
        "The correct(surfaceSymmTensorField&) function is not implemented\n"
        " for the " + type() + " mechanical law"
    );
}


Foam::scalar Foam::mechanicalLaw::residual()
{
    // Default to zero; this can be overwritten by any derived mechanical law
    return 0.0;
}


Foam::scalar Foam::mechanicalLaw::newDeltaT()
{
    // Default to a large number
    return mesh_.time().endTime().value();
}


// ************************************************************************* //
