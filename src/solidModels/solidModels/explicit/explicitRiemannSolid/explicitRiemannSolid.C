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
#include "fvcPointAverage.H"
#include "fvcInterpolate.H"
#include "globalPolyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitRiemannSolid, 0);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

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
            fvc::div((xf_ ^ tractionC_)*mesh().magSf())
        );
        am_.AMconservation(x_, rhoU_, rhoURHS, rhsRhoUAM, stage);
    }
    // if (useStabilisation_)
    // {
    //     rhoURHS +=
    //         stabilisation().stabilisation
    //         (
    //             U(),
    //             fvc::grad(U())(),
    //             (deltaT*impKf_)()
    //         );
    // }
    // if (useBulkViscosity_)
    // {
    //     rhoURHS +=
    //         fvc::div
    //         (
    //             this->mesh().Sf()*energies_.viscousPressure
    //             (
    //                 this->rho(),
    //                 fvc::interpolate(sWaveSpeed_)(),
    //                 this->gradD()
    //             )
    //         );
    // }

    surfaceScalarField rhof(fvc::interpolate(rho_));

    // Update coordinates
    x_ += rhoU_/rho_*deltaT;
    xf_ += rhoUC_/rhof*deltaT;

    if (filter_)
    {
        xN_ += pointRhoU_/mechanical().volToPoint().interpolate(rho_)*deltaT;
    }
    else
    {
        pointRhoU_ == mechanical().volToPoint().interpolate(rhoU_);
    }

     // Update linear momentum
    rhoU_ += deltaT*rhoURHS;

    // Update deformation gradient tensor
    F_ += deltaT*fvc::surfaceIntegrate(rhoUC_/rhof*mesh().Sf());
    if (pTouch_)
    {
        F_ += 0.5*(fvc::surfaceIntegrate(xf_*mesh().Sf()) - F_);
    }
    F_.correctBoundaryConditions();

    // Calculate primitive variables
    decode();
}


void explicitRiemannSolid::updateFluxes()
{
    update();

    // Surface normals
    const surfaceVectorField& N = mech_.N();
    const surfaceVectorField& n = mech_.n();

    // Reconstruction of Piola tensor
    autoPtr<ReconstructionScheme<tensor>> PLimiter
    (
        ReconstructionScheme<tensor>::New(P_, "P")
    );
    surfaceVectorField tractionOwn(PLimiter->interpolateOwn() & N);
    surfaceVectorField tractionNei(PLimiter->interpolateNei() & N);

    // Momentum
    autoPtr<ReconstructionScheme<vector>> rhoULimiter
    (
        ReconstructionScheme<vector>::New(rhoU_, "U")
    );
    surfaceVectorField rhoUOwn(rhoULimiter->interpolateOwn());
    surfaceVectorField rhoUNei(rhoULimiter->interpolateNei());

    const surfaceTensorField& stabRhoU(mech_.stabRhoU());
    const surfaceTensorField& stabTraction(mech_.stabTraction());

    // Acoustic Riemann solver
    tractionC_ =
    (
        0.5*(tractionOwn + tractionNei + (stabRhoU & (rhoUNei - rhoUOwn)))
    );
    rhoUC_ =
    (
        0.5*(rhoUOwn + rhoUNei + (stabTraction & (tractionNei - tractionOwn)))
    );

    volVectorField::Boundary& brhoU(rhoU_.boundaryFieldRef());
    surfaceVectorField::Boundary& brhoUC(rhoUC_.boundaryFieldRef());
    pointVectorField::Boundary& bpointRhoU(pointRhoU_.boundaryFieldRef());
    surfaceVectorField::Boundary& btractionC(tractionC_.boundaryFieldRef());

    pointScalarField pointRho(this->mechanical().volToPoint().interpolate(rho_));

    forAll(btractionC, patchi)
    {
        const polyPatch& p = mesh().boundaryMesh()[patchi];
        const fvPatchField<vector>& pDD(DD_.boundaryField()[patchi]);
        const vectorField& pn(n.boundaryField()[patchi]);

        if (isA<tractionBase>(pDD))
        {
            const tractionBase& tb = dynamicCast<const tractionBase>(pDD);
            vectorField tp((tb.traction() - pn*tb.pressure()));

            brhoUC[patchi] =
                rhoUOwn.boundaryField()[patchi]
              + (
                    stabTraction.boundaryField()[patchi]
                  & (tp - tractionOwn.boundaryField()[patchi])
                );
            btractionC[patchi] = tp;
        }
        else if (pDD.fixesValue())
        {
            brhoUC[patchi] ==
                rho_.boundaryField()[patchi]
               *pDD/this->mesh().time().deltaTValue();
            bpointRhoU[patchi] ==
                pointRho.boundaryField()[patchi].patchInternalField()
               *pointDD_.boundaryField()[patchi].patchInternalField()
               /this->mesh().time().deltaTValue();

            btractionC[patchi] ==
                tractionOwn.boundaryField()[patchi]
              + (
                    stabRhoU.boundaryField()[patchi]
                  & (
                        brhoUC[patchi]
                      - rhoUOwn.boundaryField()[patchi]
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
                    rhoUOwn.boundaryField()[patchi]
                  - tractionOwn.boundaryField()[patchi]
                   /sWaveSpeed_.boundaryField()[patchi]
                );
            btractionC[patchi] =
                (pn*pn)
              & (
                    tractionOwn.boundaryField()[patchi]
                  - pWaveSpeed_.boundaryField()[patchi]
                   *rhoUOwn.boundaryField()[patchi]
                );
        }
        else if (!pDD.coupled())
        {
            brhoUC[patchi] = rhoUOwn.boundaryField()[patchi];
            btractionC[patchi] = tractionOwn.boundaryField()[patchi];
        }
        brhoU[patchi] == brhoUC[patchi];
    }

    // Filter linear momentum
    if (filter_)
    {
        volVectorField rhoUAvg(fvc::surfVolInterpolate(rhoUC_));
        volTensorField rhoUGradLocal
        (
            gradSchemes_.localGradient(rhoUAvg, rhoUC_, pointRhoU_)
        );
        fvc::volPointInterpolate(rhoUAvg, rhoUGradLocal, pointRhoU_, true);
        pointRhoU_.correctBoundaryConditions();

        rhoUC_ = fvc::average(pointRhoU_);
    }
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
        solidModelDict().lookupOrDefault("conserveAngularMomentum", true)
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
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.C() + D_
    ),
    xf_
    (
        IOobject
        (
            "x",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(x_)
    ),
    xN_
    (
        IOobject
        (
            "xN",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->pointD()
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
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", rhoU_.dimensions(), Zero),
        pointDBoundaryTypes(DD_)
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
    relaxation_(this->solidModelDict().optionalSubDict("relaxation")),
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

    useStabilisation_(solidModelDict().lookupOrDefault("useStabilisation", false)),
    useBulkViscosity_(solidModelDict().lookupOrDefault("useBulkViscosity", false)),

    energies_(mesh, solidModelDict()),

    pTouch_(solidModelDict().lookupOrDefault("pTouch", false)),
    filter_(solidModelDict().lookupOrDefault("filter", false)),

    impK_("impK", mechanical().impK()),
    impKf_("impKf", mechanical().impKf()),
    curIndex_(-1)
{
    if (!xN_.headerOk())
    {
        xN_.primitiveFieldRef() += mesh.points();
        xN_.correctBoundaryConditions();
    }
    DDisRequired(type);
    if (!useStabilisation_)
    {
        stabilisation().setMethods(momentumStabilisation::NONE, dictionary());
    }
    updateFluxes();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitRiemannSolid::decode()
{
    // Enforce any cell displacements
    if (this->setCellDisps().cellIDs().size())
    {
        const scalar deltaT = this->mesh().time().deltaTValue();
        vectorField& xI = x_.primitiveFieldRef();
        vectorField& xNI = xN_.primitiveFieldRef();
        vectorField& rhoUI = rhoU_;
        vectorField& pointRhoUI = pointRhoU_;
        const vectorField& xOldI = x_.oldTime();
        const scalarField& rhoI = rho_;

        const labelList& cells = this->setCellDisps().cellIDs();
        const vectorField& cellDs = this->setCellDisps().cellDisps();
        const labelListList& cellPoints = mesh().cellPoints();
        const pointField& points = mesh().points();

        forAll(cells, i)
        {
            const label celli = cells[i];
            xI[celli] = mesh().C()[celli] + cellDs[i];
            rhoUI[celli] = rhoI[celli]*(xI[celli] - xOldI[celli])/deltaT;

            const labelList& cp = cellPoints[celli];
            forAll(cp, pi)
            {
                const label pointi = cp[pi];
                pointRhoUI[pointi] = rhoUI[celli];
                xNI[pointi] = points[pointi] + cellDs[i];
            }

        }
    }
    x_.correctBoundaryConditions();
    rhoU_.correctBoundaryConditions();

    U_ == rhoU_/rho_;

    // Update displacements
    {
        volVectorField xOld(x_);

        D_ = x_ - mesh().C();
        D_.correctBoundaryConditions();

        DD_ = x_ - xOld;
        DD_.correctBoundaryConditions();
    }

    if (filter_)
    {
        pointVectorField xNOld(xN_);
        pointD_.primitiveFieldRef() = xN_.primitiveField() - mesh().points();
        pointD_.correctBoundaryConditions();

        const pointConstraints& pcs = pointConstraints::New(pointD_.mesh());
        pcs.constrainDisplacement(pointD_, true);

        pointDD_ == xN_ - xNOld;
    }
    else
    {
        this->mechanical().volToPoint().interpolate(D_, pointD_);
        pointVectorField xNOld(xN_);
        xN_.primitiveFieldRef() = mesh().points() + pointD_.primitiveField();
        xN_.correctBoundaryConditions();

        pointDD_ == xN_ - xNOld;
    }
}


void explicitRiemannSolid::update(const bool correctSigma)
{
    //- Update gradients
    mechanical().grad(D_, gradD_);
    // gradD_ = F_ - tensor::I;
    gradDD_ = gradD_ - gradD_.oldTime();

    //- Update wavespeeds
    // Update deformation quantities
    mech_.correctDeformation();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
    impK_ = mechanical().impK();
    impKf_ = mechanical().impKf();

    pWaveSpeed_ =
        sqrt(mechanical().elasticModulus()/rho_)/beta_/mech_.stretch();
    sWaveSpeed_ =
        sqrt(mechanical().shearModulus()/rho_)*beta_/mech_.stretch();

    mech_.correct(pWaveSpeed_, sWaveSpeed_);

    P_ = mechanical().P(sigma());
    // if (useBulkViscosity_)
    // {
    //     P_ +=
    //         fvc::reconstruct
    //         (
    //             energies_.viscousPressure
    //             (
    //                 rho(),
    //                 fvc::interpolate(pWaveSpeed_),
    //                 gradD()
    //             )*mech_.n()
    //         );
    // }
}


bool explicitRiemannSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;
    const dimensionedScalar& deltaT = mesh().time().deltaT();

    mesh().update();
    enforceLinear() = false;

    // do
    {
        volVectorField rhoURHS
        (
            volVectorField::New
            (
                "rhoURHS",
                mesh(),
                dimensionedVector(rhoU_.dimensions()/dimTime, Zero)
            )
        );

        // Reset fields
        if (curIndex_ != this->runTime().timeIndex())
        {
            x_ == x_.oldTime();
            xf_ == xf_.oldTime();
            xN_ == xN_.oldTime();

            rhoU_ == rhoU_.oldTime();

            F_ == F_.oldTime();
        }
        else
        {
            curIndex_ = this->runTime().timeIndex();
        }


        // Predictor
        updateFluxes();
        solveGEqns(rhoURHS, 0);

        // Corrector
        updateFluxes();
        solveGEqns(rhoURHS, 1);


        // Average old time and new time
        x_ == 0.5*(x_.oldTime() + x_);
        xf_ == 0.5*(xf_.oldTime() + xf_);
        xN_ == 0.5*(xN_.oldTime() + xN_);
        x_.correctBoundaryConditions();
        xN_.correctBoundaryConditions();

        //- Update momentum
        rhoU_ == 0.5*(rhoU_.oldTime() + rhoU_);

        // Update deformation gradient tensor
        F_ == 0.5*(F_.oldTime() + F_);

        decode();

    }// while (mesh().update());

    // Check energies
    energies_.checkEnergies
    (
        this->rho(),
        this->U(),
        this->D(),
        this->DD(),
        this->sigma(),
        this->gradD(),
        this->gradDD(),
        this->stabilisation(),
        this->g()
    );

    {
        vector linearMoementum(Zero);
        vector angularMomentum(Zero);
        scalar totalV(0.0);

        const scalarField& V = mesh().V();
        forAll(V, celli)
        {
            linearMoementum += rhoU_[celli]*V[celli];
            angularMomentum += V[celli]*(x_[celli] ^ rhoU_[celli]);
            totalV += V[celli];
        }

        reduce(linearMoementum, sumOp<vector>());
        reduce(angularMomentum, sumOp<vector>());
        reduce(totalV, sumOp<scalar>());


        Info<< "Total linear momentum = " << linearMoementum/totalV << nl
            << "Total angular momentum = " << angularMomentum/totalV << endl;
    }

    // Re-read the dictionary
    relaxation_.read
    (
        this->solidModelDict().optionalSubDict("relaxation")
    );

    //- Relax, if wanted
    relaxation_.relax(this->U(), this->rho());
    relaxation_.relax(rhoU_);
    relaxation_.relax(pointRhoU_);

    return true;
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
