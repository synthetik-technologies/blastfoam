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

InClass
    Foam::mechanicalLaw

\*---------------------------------------------------------------------------*/

#include "mechanicalLaw.H"
#include "volFields.H"
#include "fvc.H"
#include "IOdictionary.H"
#include "lookupSolidModel.H"
#include "solidModel.H"
#include "fvm.H"
#include "fvc.H"
#include "zeroGradientFvPatchFields.H"

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
        FatalErrorIn("void " + type() + "::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_.set
    (
        new volTensorField
        (
            IOobject
            (
                "F+" + type(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeFf() const
{
    if (FfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_.set
    (
        new surfaceTensorField
        (
            IOobject
            (
                "Ff_" + type(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}

void Foam::mechanicalLaw::makeRelF() const
{
    if (relFPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeRelF()")
            << "pointer already set" << abort(FatalError);
    }

    relFPtr_.set
    (
        new volTensorField
        (
            IOobject
            (
                "relF_" + type(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        )
    );
}


void Foam::mechanicalLaw::makeRelFf() const
{
    if (relFfPtr_.valid())
    {
        FatalErrorIn("void " + type() + "::makeRelFf()")
            << "pointer already set" << abort(FatalError);
    }

    relFfPtr_.set
    (
        new surfaceTensorField
        (
            IOobject
            (
                "relFf_" + type(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
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
    if (useSolidDeformation_)
    {
        return mesh_.lookupObject<volTensorField>("F");
    }
    if (FPtr_.empty())
    {
        makeF();
    }

    return FPtr_();
}


const Foam::surfaceTensorField& Foam::mechanicalLaw::Ff() const
{
    if (useSolidDeformation_)
    {
        return mesh_.lookupObject<surfaceTensorField>("Ff");
    }

    if (FfPtr_.empty())
    {
        makeFf();
    }

    return FfPtr_();
}


const Foam::volTensorField& Foam::mechanicalLaw::relF() const
{
    if (relFPtr_.empty())
    {
        makeRelF();
    }

    return relFPtr_();
}


Foam::volTensorField& Foam::mechanicalLaw::relF()
{
    if (relFPtr_.empty())
    {
        makeRelF();
    }

    return relFPtr_();
}


const Foam::surfaceTensorField& Foam::mechanicalLaw::relFf() const
{
    if (relFfPtr_.empty())
    {
        makeRelFf();
    }

    return relFfPtr_();
}


Foam::surfaceTensorField& Foam::mechanicalLaw::relFf()
{
    if (relFfPtr_.empty())
    {
        makeRelFf();
    }

    return relFfPtr_();
}


void Foam::mechanicalLaw::makeJ() const
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


const Foam::volScalarField& Foam::mechanicalLaw::J() const
{
    if (useSolidDeformation_)
    {
        return mesh_.lookupObject<volScalarField>("J");
    }

    if (JPtr_.empty())
    {
        makeJ();
    }

    return JPtr_();
}


void Foam::mechanicalLaw::makeJf() const
{
    if (JfPtr_.valid())
    {
        FatalErrorInFunction
            << "pointer already set" << abort(FatalError);
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


const Foam::surfaceScalarField&
Foam::mechanicalLaw::Jf() const
{
    if (useSolidDeformation_)
    {
        return mesh_.lookupObject<surfaceScalarField>("Jf");
    }

    if (JfPtr_.empty())
    {
        makeJf();
    }

    return JfPtr_();
}

Foam::volScalarField& Foam::mechanicalLaw::sigmaHyd()
{
    if (sigmaHydPtr_.empty())
    {
        sigmaHydPtr_.set
        (
            new volScalarField
            (
                IOobject
                (
                    "sigmaHyd",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedScalar("zero", dimPressure, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    return sigmaHydPtr_();
}


Foam::volVectorField& Foam::mechanicalLaw::gradSigmaHyd()
{
    if (gradSigmaHydPtr_.empty())
    {
        gradSigmaHydPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    "grad(sigmaHyd)",
                    mesh().time().timeName(),
                    mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh(),
                dimensionedVector("zero", dimPressure/dimLength, Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );
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
    if (useSolidDeformation_)
    {
        const volTensorField& F =
            mesh_.lookupObject<volTensorField>("F");

        // Check if the mathematical model is in total or updated Lagrangian form
        if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            FatalErrorInFunction
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
            return false;
        }
        else
        {
            // Update the relative deformation gradient: not needed
            relF() = F & inv(F.oldTime());
            return false;
        }
    }

    if (!FPtr_.valid())
    {
        makeF();
    }
    volTensorField& F = FPtr_();

    if (!JPtr_.valid())
    {
        makeJ();
    }
    volScalarField& J = JPtr_();

    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate the relative deformation gradient
        relF() = I + gradDD.T();

        // Update the total deformation gradient
        F = relF() & F.oldTime();

        //- Calculate Jacobian
        J = det(F);

        if (enforceLinear())
        {
            WarningIn
            (
                "void Foam::MooneyRivlinThreeParametersElastic::"
                "correct(volSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime()
              + 2.0*mu*symm(gradDD) + (K - (2.0/3.0)*mu)*tr(gradDD)*I;

            return true;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const volTensorField& gradDD =
                mesh().lookupObject<volTensorField>("grad(DD)");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            F = F.oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            relF() = F & inv(F.oldTime());

            //- Calculate Jacobian
            J = det(F);

            if (enforceLinear())
            {
                WarningIn
                (
                    "void " + type() + "::correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                return true;
            }
        }
        else
        {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            // Update the total deformation gradient
            F = I + gradD.T();

            // Update the relative deformation gradient: not needed
            relF() = F & inv(F.oldTime());

            //- Calculate Jacobian
            J = det(F);

            if (enforceLinear())
            {
                WarningIn
                (
                    "void " + type() + "::correct(volSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;

                return true;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "::correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // linearised elasticity was not enforced
    return false;
}


bool Foam::mechanicalLaw::updateFf
(
    surfaceSymmTensorField& sigma,
    const dimensionedScalar& mu,
    const dimensionedScalar& K
)
{
    if (useSolidDeformation_)
    {
        const surfaceTensorField& Ff =
            mesh_.lookupObject<surfaceTensorField>("Ff");

        // Check if the mathematical model is in total or updated Lagrangian form
        if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            FatalErrorInFunction
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
            return false;
        }
        else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            if (incremental())
            {
                // Update the relative deformation gradient: not needed
                relFf() = Ff & inv(Ff.oldTime());

                if (enforceLinear())
                {
                    WarningInFunction
                        << "Material linearity enforced for stability!"
                        << endl;

                    // Lookup gradient of displacement increment
                    const surfaceTensorField& gradDD =
                        mesh().lookupObject<surfaceTensorField>("grad(DD)");

                    // Calculate stress using Hooke's law
                    sigma =
                        sigma.oldTime()
                    + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                    return true;
                }
            }
            else
            {
                // Update the relative deformation gradient: not needed
                relFf() = Ff & inv(Ff.oldTime());

                if (enforceLinear())
                {
                    WarningInFunction
                        << "Material linearity enforced for stability!"
                        << endl;

                    // Lookup gradient of displacement
                    const surfaceTensorField& gradD =
                        mesh().lookupObject<surfaceTensorField>("grad(D)");

                    // Calculate stress using Hooke's law
                    sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;

                    return true;
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unknown nonLinGeom type: " << nonLinGeom()
                << abort(FatalError);
            return false;
        }
    }

    if (!FfPtr_.valid())
    {
        makeFf();
    }
    surfaceTensorField& Ff = FfPtr_();

    if (!JfPtr_.valid())
    {
        makeJf();
    }
    surfaceScalarField& Jf = JfPtr_();

    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(surfaceSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient: not needed
        relFf() = I + gradDD.T();

        // Update the total deformation gradient
        Ff = relFf() & Ff.oldTime();

        //- Calculate Jacobian
        Jf = det(Ff);

        if (enforceLinear())
        {
            WarningIn
            (
                "void " + type() + "::correct(surfaceSymmTensorField& sigma)"
            )   << "Material linearity enforced for stability!" << endl;

            // Calculate stress using Hooke's law
            sigma =
                sigma.oldTime() + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

            return true;
        }
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        if (incremental())
        {
            // Lookup gradient of displacement increment
            const surfaceTensorField& gradDD =
                mesh().lookupObject<surfaceTensorField>("grad(DD)f");

            // Update the total deformation gradient
            // Note: grad is wrt reference configuration
            Ff = Ff.oldTime() + gradDD.T();

            // Update the relative deformation gradient: not needed
            relFf() = Ff & inv(Ff.oldTime());

            //- Calculate Jacobian
            Jf = det(Ff);

            if (enforceLinear())
            {
                WarningIn
                (
                    "void " + type()
                  + "::correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma =
                    sigma.oldTime()
                  + 2.0*mu*dev(symm(gradDD)) + K*tr(gradDD)*I;

                return true;
            }
        }
        else
        {
            // Lookup gradient of displacement
            const surfaceTensorField& gradD =
                mesh().lookupObject<surfaceTensorField>("grad(D)f");

            // Update the total deformation gradient
            Ff = I + gradD.T();

            // Update the relative deformation gradient: not needed
            relFf() = Ff & inv(Ff.oldTime());

            //- Calculate Jacobian
            Jf = det(Ff);

            if (enforceLinear())
            {
                WarningIn
                (
                    "void " + type()
                  + "::correct(surfaceSymmTensorField& sigma)"
                )   << "Material linearity enforced for stability!" << endl;

                // Calculate stress using Hooke's law
                sigma = 2.0*mu*dev(symm(gradD)) + K*tr(gradD)*I;



                return true;
            }
        }
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "::correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // linearised elasticity was not enforced
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
        sigmaHyd().storePrevIter();

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
                "void " + type() + "updateSigmaHyd(...)\n"
            )   << "Cannot find the DEqnA or DDEqnA field: this should be "
                << "stored in the solidModel" << abort(FatalError);
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
        sigmaHyd().relax();
    }
    else
    {
        // Explicitly calculate hydrostatic stress
        // We use 1.0* to overwritting the field IOobject attributes e.g. its
        // name and writeOpt
        sigmaHyd() = 1.0*sigmaHydExplicit;
    }

    // Update the gradient
    gradSigmaHyd() = fvc::grad(sigmaHyd());
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
    rho_(mesh.lookupObject<volScalarField>("rho")),
    FPtr_(),
    FfPtr_(),
    relFPtr_(),
    relFfPtr_(),
    sigmaHydPtr_(),
    gradSigmaHydPtr_(),
    useSolidDeformation_(false),
    solvePressureEqn_
    (
        dict.lookupOrDefault<Switch>("solvePressureEqn", false)
    ),
    pressureSmoothingScaleFactor_
    (
        dict.lookupOrDefault<scalar>("pressureSmoothingScaleFactor", 100.0)
    )
{
    // Set the base mesh region name
    // For an FSI case, the region will be called solid, else it will be called
    // region0.
//     if (mesh.time().foundObject<fvMesh>("solid"))
//     {
//         baseMeshRegionName_ = "solid";
//     }
//     else if (mesh.time().foundObject<fvMesh>("region0"))
//     {
//         baseMeshRegionName_ = "region0";
//     }
//     else
//     {
//         FatalErrorIn
//         (
//             "Foam::mechanicalLaw::mechanicalLaw\n"
//             "(\n"
//             "    const word& name,\n"
//             "    const fvMesh& mesh,\n"
//             "    const dictionary& dict\n"
//             ")"
//         ) << "solid region name not found" << abort(FatalError);
//     }

    if (solvePressureEqn_)
    {
        Info<< "    Laplacian equation will be solved for pressure" << nl
            << "    pressureSmoothingScaleFactor: "
            << pressureSmoothingScaleFactor_
            << endl;
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
