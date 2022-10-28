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

#include "explicitRiemannSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "fvcCellReduce.H"
#include "labelVector.H"
#include "addToRunTimeSelectionTable.H"
#include "primitivePatchInterpolation.H"

#include "meshSizeObject.H"

#include "wedgePolyPatch.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "solidTractionFvPatchVectorField.H"
#include "ReconstructionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitRiemannSolid, 0);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void explicitRiemannSolid::updateStress()
{
    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
    impKf_ = mechanical().impKf();

    const volScalarField& J = mech_.J();
    const volTensorField& invF = mech_.invF();
    const volSymmTensorField& sigma = this->sigma();
    forAll(P_, celli)
    {
        P_[celli] = J[celli]*(invF[celli] & sigma[celli]);
    }
    volTensorField::Boundary& bP = P_.boundaryFieldRef();
    forAll(bP, patchi)
    {
        fvPatchTensorField& pP = bP[patchi];
        const fvPatchSymmTensorField& psigma = sigma.boundaryField()[patchi];
        const fvPatchScalarField& pJ = J.boundaryField()[patchi];
        const fvPatchTensorField& pinvF = invF.boundaryField()[patchi];
        forAll(pP, facei)
        {
            pP[facei] = pJ[facei]*(pinvF[facei] & psigma[facei]);
        }
    }
}


void explicitRiemannSolid::solveGEqns
(
    volVectorField& rhoURHS,
    const label stage
)
{
    dimensionedScalar deltaT(mesh().time().deltaT());

    // Compute right hand sides
    rhoURHS = fvc::surfaceIntegrate(tractionC_*mesh().magSf());
    if (angularMomentumConservation_)
    {
        volVectorField rhsRhoUAM
        (
            fvc::div(((Df_ + mesh().Cf()) ^ tractionC_)*mesh().magSf())
        );
        am_.AMconservation(x_, rhoU_, rhoURHS, rhsRhoUAM, stage);
    }
    if (useBulkViscosity_)
    {
        rhoURHS +=
            fvc::div
            (
                mesh().Sf()*energies_.viscousPressure
                (
                    rho(), fvc::interpolate(pWaveSpeed_), gradD()
                )
            );
    }
    if (useStabilisation_)
    {
        rhoURHS +=
            stabilisation().stabilisation
            (
                U(),
                fvc::grad(U())(),
                (deltaT*impKf_)()
            );
    }

    if (mesh().relaxField(D_.name()))
    {
        scalar fac = mesh().fieldRelaxationFactor(D_.name());
        rhoURHS *= fac;
    }

    surfaceScalarField rhof(fvc::interpolate(rho_));

    // Update coordinates
    {
        const volVectorField prevD(D_);
        const volTensorField prevGradD(gradD_);
        const pointVectorField prevPointD(pointD_);

        D_ += U_*deltaT;
        D_.correctBoundaryConditions();

        Df_ += deltaT*rhoUC_/rhof;
        pointD_ += deltaT*pointRhoU_/mechanical().volToPoint().interpolate(rho_);

        // mechanical().interpolate(D_, pointD_);
        mechanical().grad(D_, gradD_);

        DD_ == D_ - prevD;
        gradDD_ = gradD_ - prevGradD;
        pointDD_ = pointD_ - prevPointD;
    }

    // Face displacement
    Df_ = fvc::interpolate(D_);

    // Material positions
    x_ = mesh().C() + D_;
    x_.correctBoundaryConditions();

    // Update linear momentum
    rhoU_ += deltaT*rhoURHS;

    U_ = rhoU_/rho_;
    U_.boundaryFieldRef() == DD_.boundaryField()/deltaT.value();
    rhoU_.boundaryFieldRef() == U_.boundaryField()*rho_.boundaryField();

    // Update deformation gradient tensor
    F_ += deltaT*fvc::surfaceIntegrate(rhoUC_/rhof*mesh().Sf());
    // F_ = I + gradD_.T();

    // Update deformation quantities
    mech_.correctDeformation();
}


void explicitRiemannSolid::updateFluxes()
{
    pWaveSpeed_ =
        sqrt(mechanical().elasticModulus()/rho_)/beta_/mech_.stretch();
    sWaveSpeed_ =
        sqrt(mechanical().shearModulus()/rho_)*beta_/mech_.stretch();

    mech_.correct(pWaveSpeed_, sWaveSpeed_);

    updateStress();

    // Surface normals
    const surfaceVectorField& N = mech_.N();
    const surfaceVectorField& n = mech_.n();

    autoPtr<ReconstructionScheme<tensor>> PLimiter
    (
        ReconstructionScheme<tensor>::New(P_, "P")
    );

    // Reconstruction
    surfaceTensorField POwn
    (
        surfaceTensorField::New
        (
            "POwn",
            PLimiter->interpolateOwn()
        )
    );
    surfaceTensorField PNei
    (
        surfaceTensorField::New
        (
            "PNei",
            PLimiter->interpolateNei()
        )
    );
    tractionOwn_ = POwn & N;
    tractionNei_ = PNei & N;

    autoPtr<ReconstructionScheme<vector>> rhoULimiter
    (
        ReconstructionScheme<vector>::New(rhoU_, "U")
    );
    rhoUOwn_ = rhoULimiter->interpolateOwn();
    rhoUNei_ = rhoULimiter->interpolateNei();

    const surfaceTensorField& stabRhoU(mech_.stabRhoU());
    const surfaceTensorField& stabTraction(mech_.stabTraction());

    // Acoustic Riemann solver
    tractionC_ =
        0.5*(tractionOwn_ + tractionNei_)
      + 0.5*(stabRhoU & (rhoUNei_ - rhoUOwn_));
    rhoUC_ =
        0.5*(rhoUOwn_ + rhoUNei_)
      + 0.5*(stabTraction & (tractionNei_ - tractionOwn_));

    surfaceVectorField::Boundary& brhoUC(rhoUC_.boundaryFieldRef());
    surfaceVectorField::Boundary& btractionC(tractionC_.boundaryFieldRef());
    forAll(btractionC, patchi)
    {
        const polyPatch& p = mesh().boundaryMesh()[patchi];
        const fvPatch& patch = mesh().boundary()[patchi];
        const fvPatchField<vector>& pD(D_.boundaryField()[patchi]);
        const fvPatchField<vector>& pRhoU(rhoU_.boundaryField()[patchi]);
        const vectorField& pn(n.boundaryField()[patchi]);

        // Riemann solver for inter-processor boundaries
        if (isA<solidTractionFvPatchVectorField>(pD))
        {
            const solidTractionFvPatchVectorField& stD =
                dynamicCast<const solidTractionFvPatchVectorField>(pD);
            vectorField tp
            (
                (stD.traction() - pn*stD.pressure())
            );

            brhoUC[patchi] =
                rhoUOwn_.boundaryField()[patchi]
              + (
                    stabTraction.boundaryField()[patchi]
                  & (tp - tractionOwn_.boundaryField()[patchi])
                );
            btractionC[patchi] = tp;
        }
        else if (pD.fixesValue() || U_.boundaryField()[patchi].fixesValue())
        {
            brhoUC[patchi] = pRhoU[patchi];

            btractionC[patchi] =
                tractionOwn_.boundaryField()[patchi]
              + (
                    stabRhoU.boundaryField()[patchi]
                  & (
                        brhoUC[patchi]
                      - rhoUOwn_.boundaryField()[patchi]
                    )
                );
        }
        else if
        (
            isA<symmetryPolyPatch>(p)
         || isA<symmetryPlanePolyPatch>(p)
        )
        {
            brhoUC[patchi] =
                (tensor::I - pn*pn)
              & (
                    rhoUOwn_.boundaryField()[patchi]
                  - tractionOwn_.boundaryField()[patchi]
                   /sWaveSpeed_.boundaryField()[patchi]
                );
            btractionC[patchi] =
                (pn*pn)
              & (
                    tractionOwn_.boundaryField()[patchi]
                  - pWaveSpeed_.boundaryField()[patchi]
                   *rhoUOwn_.boundaryField()[patchi]
                );
        }
    }

    // Average linear momentum
    volVectorField rhoUAvg(fvc::average(rhoUC_));
    pointRhoU_ = volPointInterpolation::New(mesh()).interpolate(rhoUAvg);

    // rhoUC_.ref() = interpSchemes_.pointToSurface(pointRhoU_)()();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitRiemannSolid::explicitRiemannSolid
(
    const word& type,
    dynamicFvMesh& mesh,
    const nonLinearGeometry::nonLinearType& nonLinear
)
:
    solidModel(type, mesh, nonLinear, incremental()),
    F_
    (
        IOobject
        (
            "F",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("I", dimless, I)
    ),
    ops_(mesh),
    beta_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "incompressiblilityCoefficient",
            1.0
        )
     ),
    angularMomentumConservation_
    (
        solidModelDict().lookupOrDefault("angularMomentumConservation", true)
    ),
    D_(this->D()),
    DD_(this->DD()),
    gradD_(this->gradD()),
    gradDD_(this->gradDD()),
    x_
    (
        IOobject
        (
            "x",
            mesh.time().timeName(),
            mesh
        ),
        mesh.C() + D_
    ),
    Df_
    (
        IOobject
        (
            "Df",
            mesh.time().timeName(),
            mesh
        ),
        fvc::interpolate(D_)
    ),
    pointD_(this->pointD()),
    pointDD_(this->pointDD()),
    rho_(this->rho()),
    U_(this->U()),
    rhoU_
    (
        IOobject
        (
            "rhoU",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        rho_*U_
    ),
    rhoUC_
    (
        IOobject
        (
            "rhoUC",
            mesh.time().timeName(),
            mesh
        ),
        fvc::interpolate(rho_*U_)
    ),
    pointRhoU_
    (
        IOobject
        (
            "pointRhoU",
            mesh.time().timeName(),
            mesh
        ),
        pMesh(),
        dimensionedVector("0", rhoU_.dimensions(), Zero)
        // pointDBoundaryTypes(D_)
    ),
    tractionC_
    (
        IOobject
        (
            "tractionC",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("0", dimensionSet(1, -1, -2, 0, 0, 0, 0), Zero)
    ),
    mech_(F_, ops_),
    interpSchemes_(mesh),
    gradSchemes_(DD_),
    am_(mesh, *this),
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
        mech_.J()*(sigma() & ops_.invT(F_))
    ),
    pWaveSpeed_
    (
        IOobject
        (
            "pWaveSpeed",
            mesh.time().timeName(),
            mesh
        ),
        sqrt(mechanical().elasticModulus()/rho_)/beta_/mech_.stretch()
    ),
    sWaveSpeed_
    (
        IOobject
        (
            "sWaveSpeed",
            mesh.time().timeName(),
            mesh
        ),
        sqrt(mechanical().shearModulus()/rho_)*beta_/mech_.stretch()
    ),
    rhoUOwn_("rhoUOwn", rhoUC_),
    rhoUNei_("rhoUNei", rhoUC_),
    tractionOwn_("tractionOwn", tractionC_),
    tractionNei_("tractionNei", tractionC_),

    useStabilisation_(solidModelDict().lookupOrDefault("useStabilisation", true)),
    useBulkViscosity_(solidModelDict().lookupOrDefault("useBulkViscosity", true)),

    energies_(mesh, solidModelDict()),
    impKf_(mechanical().impKf())
{
    DisRequired(type);
    if (!useStabilisation_)
    {
        stabilisation().setMethods(momentumStabilisation::NONE, dictionary());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool explicitRiemannSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    label iter = 0;
    do
    {
        if (iter++ > 0)
        {
            // Momentum
            rhoU_ = rhoU_.oldTime();

            // Displacements
            D_ = D_.oldTime();
            Df_ = Df_.oldTime();
            pointD_ = pointD_.oldTime();

            // Deformation
            gradD_ = gradD_.oldTime();
            F_ = F_.oldTime();
        }

        volVectorField rhoURHS
        (
            volVectorField::New
            (
                "rhoURHS",
                mesh(),
                dimensionedVector(rhoU_.dimensions()/dimTime, Zero)
            )
        );

        // Predictor
        updateFluxes();
        solveGEqns(rhoURHS, 0);

        // Corrector
        updateFluxes();
        solveGEqns(rhoURHS, 1);

        // Update coordinates
        rhoU_ = 0.5*(rhoU_.oldTime() + rhoU_);
        U_ = rhoU_/rho_;
        U_.correctBoundaryConditions();

        D_ = 0.5*(D_.oldTime() + D_);
        Df_ = 0.5*(Df_.oldTime() + Df_);
        pointD_ = 0.5*(pointD_.oldTime() + pointD_);

        D_.correctBoundaryConditions();
        DD_ = D_ - D_.oldTime();
        pointDD_ = pointD_ - pointD_.oldTime();

        // Update coordinates
        x_ = mesh().C() + D_;

        U_.boundaryFieldRef() = DD_.boundaryField()/deltaT.value();
        rhoU_.boundaryFieldRef() = U_.boundaryField()*rho_.boundaryField();

        // Update gradient of displacement increment
        mechanical().grad(D_, gradD_);

        // Update the gradient of total displacement
        gradDD_ = gradD_ - gradD_.oldTime();

        // Update deformation gradient tensor
        F_ = 0.5*(F_.oldTime() + F_);

        // Update deformation quantities
        mech_.correctDeformation(true);

    } while (mesh().update());

    return true;
}


tmp<vectorField> explicitRiemannSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField& impK(impKf_.boundaryField()[patch.index()]);

    // Patch gradient
    const tensorField& pGradD = gradD_.boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n(patch.nf());

    // Patch unit normals (deformed configuration)
    const vectorField& nCurrent(mech_.n().boundaryField()[patchID]);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )/impK
        )
    );
}


void explicitRiemannSolid::setDeltaT(Time& runTime)
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value
    const scalar requiredDeltaT =
        1.0
       /gMax
        (
            (
                mesh().surfaceInterpolation::deltaCoeffs()
               *fvc::interpolate(pWaveSpeed_)
            )()
        );

    // Lookup the desired Courant number
    const scalar maxCo =
        runTime.controlDict().lookupOrDefault<scalar>("maxCo", 0.7071);

    const scalar newDeltaT = maxCo*requiredDeltaT;

    Info<< "maxCo = " << maxCo << nl
        << "deltaT = " << newDeltaT << nl << endl;

    runTime.setDeltaT(newDeltaT);
}


scalar explicitRiemannSolid::CoNum() const
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value

    surfaceScalarField amaxSf(fvc::interpolate(pWaveSpeed_)*mesh().magSf());

    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgePolyPatch>(mesh().boundaryMesh()[patchi]))
        {
            amaxSf.boundaryFieldRef() = Zero;
        }
    }
    scalarField sumAmaxSf
    (
        fvc::surfaceSum(amaxSf)().primitiveField()
    );
    return 0.5*gMax(sumAmaxSf/mesh().V().field())*mesh().time().deltaTValue();
}


scalar explicitRiemannSolid::maxCoNum() const
{
    return
        mesh().time().controlDict().lookupOrDefault<scalar>
        (
            "maxCo",
            0.7071
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
