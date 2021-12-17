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

#include "explicitTotalLagSolid.H"
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitTotalLagSolid, 0);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void explicitTotalLagSolid::updateStress()
{
    // Update gradient of displacement
    mechanical().grad(D_, gradD());

    // Update gradient of displacement increment
    mechanical().grad(DD(), gradDD());

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
    impK_ = mechanical().impK();

    P_ = mech_.J()*(sigma() & ops_.invT(F_));
}


void explicitTotalLagSolid::solveGEqns
(
    volVectorField& rhoURHS,
    volVectorField& rhoURHS1,
    const label stage
)
{
    dimensionedScalar deltaT(mesh().time().deltaT());

    // Compute right hand sides
    rhoURHS = fvc::surfaceIntegrate(tractionC_*mesh().magSf());
    if (JSTScaleFactor_ > 0)
    {
        dimensionedScalar deltaT0(mesh().time().deltaT0());
        rhoURHS -=
            JSTScaleFactor_*fvc::laplacian
            (
                mesh().magSf(),
                fvc::laplacian
                (
                    deltaT*mechanical().impKf(),
                    U(),
                    "laplacian(DU,U)"
                ),
                "laplacian(DU,U)"
            );
    }
    if (angularMomentumConservation_)
    {
        volVectorField rhsRhoUAM
        (
            fvc::div(((Df_ + mesh().Cf()) ^ tractionC_)*mesh().magSf())
        );
        am_.AMconservation(rhoURHS, rhoURHS1, rhsRhoUAM, stage);
    }
    vector validD((mesh().geometricD() + labelVector::one)/2);
    rhoURHS = cmptMultiply(validD, rhoURHS);

    // Update coordinates
    volVectorField prevD(D_);
    D_ += deltaT*rhoU_/rho();
    D_.correctBoundaryConditions();

    Df_ += deltaT*rhoUC_/fvc::interpolate(rho());

    pointVectorField prevPointD(pointD_);
    pointScalarField pointRho
    (
        volPointInterpolation::New(mesh()).interpolate(rho())
    );

    pointD_ += deltaT*pointRhoU_/pointRho;
    const pointConstraints& pcs = pointConstraints::New(pointD_.mesh());
    pcs.constrainDisplacement(pointD_, false);
    pointRhoU_ = pointRho*(pointD_ - prevPointD)/deltaT;

    // Update increment of displacement
    x_ = mesh().C() + D_;
    DD() = D_ - prevD;
    pointDD() = pointD_ - prevPointD;

    // Update linear momentum
    rhoU_ += deltaT*rhoURHS;
    U_.ref() = rhoU_()/rho()();
    U_.boundaryFieldRef() = DD().boundaryField()/deltaT.value();
    rhoU_.boundaryFieldRef() = rho().boundaryField()*U_.boundaryField();

    // Update deformation gradient tensor
    F_ += deltaT*fvc::div(rhoUC_/fvc::interpolate(rho())*mesh().Sf());
}


void explicitTotalLagSolid::updateFluxes()
{
    mech_.correct(pWaveSpeed_, sWaveSpeed_);

    pWaveSpeed_ =
        sqrt(mechanical().elasticModulus()/rho())/beta_/mech_.stretch();
    sWaveSpeed_ =
        sqrt(mechanical().shearModulus()/rho())*beta_/mech_.stretch();

    updateStress();

    // Surface normals
    const surfaceVectorField& N = mech_.N();
    const surfaceVectorField& n = mech_.n();

    volVectorField Px(ops_.decomposeTensorX(P_));
    volVectorField Py(ops_.decomposeTensorY(P_));
    volVectorField Pz(ops_.decomposeTensorZ(P_));

    volTensorField gradRhoU(fvc::grad(rhoU_));
    volTensorField gradPx(fvc::grad(Px));
    volTensorField gradPy(fvc::grad(Py));
    volTensorField gradPz(fvc::grad(Pz));

    // Reconstruction
    surfaceVectorField rhoUOwn
    (
        surfaceVectorField::New
        (
            "rhoUOwn",
            mesh(),
            rhoU_.dimensions()
        )
    );
    surfaceVectorField rhoUNei
    (
        surfaceVectorField::New
        (
            "rhoUNei",
            mesh(),
            rhoU_.dimensions()
        )
    );
    surfaceTensorField POwn
    (
        surfaceTensorField::New
        (
            "POwn",
            mesh(),
            P_.dimensions()
        )
    );
    surfaceTensorField PNei
    (
        surfaceTensorField::New
        (
            "PNei",
            mesh(),
            P_.dimensions()
        )
    );


    gradSchemes_.reconstruct(rhoU_, gradRhoU, rhoUOwn, rhoUNei);
    gradSchemes_.reconstruct(P_, gradPx, gradPy, gradPz, POwn, PNei);
    surfaceVectorField tractionOwn("tractionOwn", POwn & N);
    surfaceVectorField tractionNei("tractionNei", PNei & N);

    const surfaceTensorField& stabRhoU(mech_.stabRhoU());
    const surfaceTensorField& stabTraction(mech_.stabTraction());

    // Acoustic Riemann solver
    tractionC_.ref() =
        0.5*(tractionOwn() + tractionNei())
      + 0.5*(stabRhoU() & (rhoUNei() - rhoUOwn()));
    rhoUC_.ref() =
        0.5*(rhoUOwn() + rhoUNei())
      + 0.5*(stabTraction() & (tractionNei() - tractionOwn()));

    D_.correctBoundaryConditions();
    rhoU_.correctBoundaryConditions();
    P_.correctBoundaryConditions();

    surfaceVectorField::Boundary& prhoUC(rhoUC_.boundaryFieldRef());
    surfaceVectorField::Boundary& ptractionC(tractionC_.boundaryFieldRef());
    pointVectorField::Boundary& ppointRhoU(pointRhoU_.boundaryFieldRef());

    pointScalarField pointRho(volPointInterpolation::New(mesh()).interpolate(rho()));
    forAll(ptractionC, patchi)
    {
        const polyPatch& p = mesh().boundaryMesh()[patchi];
        const fvPatch& patch = mesh().boundary()[patchi];
        const fvPatchField<vector>& pD(D_.boundaryField()[patchi]);
        const vectorField pn(n.boundaryField()[patchi]);
        const vectorField pN(N.boundaryField()[patchi]);

        // Riemann solver for inter-processor boundaries
        if (pD.fixesValue() && U_.boundaryField()[patchi].fixesValue())
        {
            prhoUC[patchi] =
                rho().boundaryField()[patchi]
               *DD().boundaryField()[patchi]
              /mesh().time().deltaT0().value();

            ptractionC[patchi] =
                tractionOwn.boundaryField()[patchi]
              + (
                    stabRhoU.boundaryField()[patchi]
                  & (
                        prhoUC[patchi]
                      - rhoUOwn.boundaryField()[patchi]
                    )
                );
            ppointRhoU[patchi].setInInternalField
            (
                pointRhoU_.primitiveFieldRef(),
                globalPatches()[p].interpolator().faceToPointInterpolate
                (
                    prhoUC[patchi]
                )()
            );
        }
        else if (isA<solidTractionFvPatchVectorField>(pD))
        {
            const solidTractionFvPatchVectorField& stD =
                dynamicCast<const solidTractionFvPatchVectorField>(pD);
            vectorField tp(stD.traction() - stD.pressure()*pn);

            prhoUC[patchi] =
                rhoUOwn.boundaryField()[patchi]
              + (
                    stabTraction.boundaryField()[patchi]
                  & (
                        tp - tractionOwn.boundaryField()[patchi]
                    )
                );
            ptractionC[patchi] = tp;
        }
        else if (patch.coupled())
        {
            const vectorField prhoUOwn
            (
                rhoU_.boundaryField()[patchi].patchInternalField()
            );
            const vectorField prhoUNei
            (
                rhoU_.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pPOwn
            (
                P_.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pPNei
            (
                P_.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pPxOwn
            (
                Px.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pPxNei
            (
                Px.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pPyOwn
            (
                Py.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pPyNei
            (
                Py.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pPzOwn
            (
                Pz.boundaryField()[patchi].patchInternalField()
            );
            const vectorField pPzNei
            (
                Pz.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradRhoUOwn
            (
                gradRhoU.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradRhoUNei
            (
                gradRhoU.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradPxOwn
            (
                gradPx.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradPxNei
            (
                gradPx.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradPyOwn
            (
                gradPy.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradPyNei
            (
                gradPy.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradPzOwn
            (
                gradPz.boundaryField()[patchi].patchInternalField()
            );
            const tensorField pgradPzNei
            (
                gradPz.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pdeltaOwn(patch.fvPatch::delta());
            const vectorField pdeltaNei(pdeltaOwn - patch.delta());

            const scalarField ppWaveSpeedOwn
            (
                pWaveSpeed_.boundaryField()[patchi].patchInternalField()
            );
            const scalarField ppWaveSpeedNei
            (
                pWaveSpeed_.boundaryField()[patchi].patchNeighbourField()
            );

            const scalarField psWaveSpeedOwn
            (
                sWaveSpeed_.boundaryField()[patchi].patchInternalField()
            );
            const scalarField psWaveSpeedNei
            (
                sWaveSpeed_.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(mesh().boundary()[patchi], facei)
            {
                const vector& dOwn(pdeltaOwn[facei]);
                const vector& dNei(pdeltaNei[facei]);

                const vector rhoUOwni =
                    prhoUOwn[facei] + (pgradRhoUOwn[facei] & dOwn);
                const vector rhoUNeii =
                    prhoUNei[facei] + (pgradRhoUNei[facei] & dNei);

                const vector PxOwni =
                    pPxOwn[facei] + (pgradPxOwn[facei] & dOwn);
                const vector PxNeii =
                    pPxNei[facei] + (pgradPxNei[facei] & dNei);

                const vector PyOwni =
                    pPyOwn[facei] + (pgradPyOwn[facei] & dOwn);
                const vector PyNeii =
                    pPyNei[facei] + (pgradPyNei[facei] & dNei);

                const vector PzOwni =
                    pPzOwn[facei] + (pgradPzOwn[facei] & dOwn);
                const vector PzNeii =
                    pPzNei[facei] + (pgradPzNei[facei] & dNei);

                const tensor POwni = tensor(PxOwni, PyOwni, PzOwni);
                const tensor PNeii = tensor(PxNeii, PyNeii, PzNeii);

                const scalar pWSi =
                    0.5*(ppWaveSpeedOwn[facei] + ppWaveSpeedNei[facei]);
                const scalar sWSi =
                    0.5*(psWaveSpeedOwn[facei] + psWaveSpeedNei[facei]);

                const vector& Ni = pN[facei];
                const vector& ni = pn[facei];

                const tensor stabRhoUi =
                    (pWSi*ni*ni) + (sWSi*(I - (ni*ni)));
                const tensor stabTractioni =
                    ((ni*ni)/pWSi) + ((I - (ni*ni))/sWSi);

                ptractionC[patchi][facei] =
                    0.5*((POwni + PNeii) & Ni)
                  + (0.5*stabRhoUi & (rhoUNeii - rhoUOwni));

                prhoUC[patchi][facei] =
                    0.5*(rhoUOwni + rhoUNeii)
                  + 0.5*(stabTractioni & ((PNeii - POwni) & Ni));
            }
        }
        else if
        (
            isA<symmetryPolyPatch>(p)
         || isA<symmetryPlanePolyPatch>(p)
        )
        {
            prhoUC[patchi] =
                (tensor::I - pn*pn)
              & (
                    rhoUOwn.boundaryField()[patchi]
                  - tractionOwn.boundaryField()[patchi]
                   /sWaveSpeed_.boundaryField()[patchi]
                );
            ptractionC[patchi] =
                (pn*pn)
              & (
                    tractionOwn.boundaryField()[patchi]
                  - pWaveSpeed_.boundaryField()[patchi]
                   *rhoUOwn.boundaryField()[patchi]
                );
            forAll(p, facei)
            {
                const label& faceID = p.start() + facei;
                forAll(mesh().faces()[faceID], node)
                {
                    const label& nodeID = mesh().faces()[faceID][node];
                    pointRhoU_[nodeID] =
                        (tensor::I - (pN[facei]*pN[facei]))
                      & pointRhoU_[nodeID];
                }
            }
        }
        else if (!polyPatch::constraintType(p.type()))
        {
            prhoUC[patchi] =
                rho().boundaryField()[patchi]
               *(
                   D_.boundaryField()[patchi]
                 - D_.oldTime().boundaryField()[patchi]
                )/mesh().time().deltaT0().value();

            ptractionC[patchi] =
                tractionOwn.boundaryField()[patchi]
              + (
                    stabRhoU.boundaryField()[patchi]
                  & (
                        prhoUC[patchi]
                      - rhoUOwn.boundaryField()[patchi]
                    )
                );
        }
        else
        {
            prhoUC[patchi] = rhoU_.boundaryField()[patchi];
            ptractionC[patchi] =  P_.boundaryField()[patchi] & pn;
        }
    }

    // Nodal linear momentum
    volVectorField rhoUAvg
    (
        fvc::average(rhoUC_)
//         interpSchemes_.surfaceToVol(rhoUC_, pointRhoU_)
    );
    volTensorField gradRhoUAvg
    (
//         fvc::grad(rhoUAvg)
        gradSchemes_.localGradient(rhoUAvg, rhoUC_, pointRhoU_)
    );
//     pointRhoU_ = volPointInterpolation::New(mesh()).interpolate(rhoUAvg);
    interpSchemes_.volToPoint(rhoUAvg, gradRhoUAvg, pointRhoU_);
    volPointInterpolation::New(mesh()).interpolateBoundaryField
    (
        rhoUAvg,
        pointRhoU_
    );

    // Symmetric boundary patch
    pointRhoU_.correctBoundaryConditions();

    // Constrained fluxes
    rhoUC_.ref() = interpSchemes_.pointToSurface(pointRhoU_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitTotalLagSolid::explicitTotalLagSolid
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
        solidModelDict().lookup<scalar>
        (
            "incompressiblilityCoefficient"
        )
     ),
    angularMomentumConservation_
    (
        solidModelDict().lookup("angularMomentumConservation")
    ),
    D_(this->D()),
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
    Df_
    (
        IOobject
        (
            "Df",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimLength, Zero)
    ),
    pointD_(this->pointD()),
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
        rho()*U()
    ),
    rhoUC_
    (
        IOobject
        (
            "rhoUC",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(rho()*U_)
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
        dimensionedVector("0", rhoU_.dimensions(), Zero)
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
    gradSchemes_(mesh),
    am_(mesh, *this),
    P_
    (
        IOobject
        (
            "P",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqrt(mechanical().elasticModulus()/rho())/beta_/mech_.stretch()
    ),
    sWaveSpeed_
    (
        IOobject
        (
            "sWaveSpeed",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqrt(mechanical().shearModulus()/rho())*beta_/mech_.stretch()
    ),
    JSTScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("JSTScaleFactor", 0.0)
    ),
    energies_(mesh, solidModelDict()),
    impK_(mechanical().impK())
{
    mechanical().setUseSolidDeformation();

    DisRequired();
    this->displacementFromVelocity(D_, DD(), 1);

    // Nodal linear momentum
    volVectorField rhoUAvg
    (
        interpSchemes_.surfaceToVol(rhoUC_, pointRhoU_)
    );
    volTensorField gradRhoUAvg
    (
        gradSchemes_.localGradient(rhoUAvg, rhoUC_, pointRhoU_)
    );
    interpSchemes_.volToPoint(rhoUAvg, gradRhoUAvg, pointRhoU_);

    // Symmetric boundary patch
    pointRhoU_.correctBoundaryConditions();

    // Constrained fluxes
    rhoUC_ = interpSchemes_.pointToSurface(pointRhoU_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitTotalLagSolid::setDeltaT(Time& runTime)
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
        1.0/
        gMax
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


bool explicitTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    rhoU_.oldTime();
    F_.oldTime();
    D_.oldTime();
    Df_.oldTime();
    pointD_.oldTime();

    volVectorField rhoURHS
    (
        volVectorField::New
        (
            "rhoURHS",
            mesh(),
            dimensionedVector(rhoU_.dimensions()/dimTime, Zero)
        )
    );
    volVectorField rhoURHS1("rhoURHS1", rhoURHS);

    updateFluxes();
    solveGEqns(rhoURHS, rhoURHS1, 0);

    updateFluxes();
    solveGEqns(rhoURHS, rhoURHS1, 1);

    rhoU_ = 0.5*(rhoU_.oldTime() + rhoU_);
    F_ = 0.5*(F_.oldTime() + F_);
    D_ = 0.5*(D_.oldTime() + D_);
    Df_ = 0.5*(Df_.oldTime() + Df_);
    pointD_ = 0.5*(pointD_.oldTime() + pointD_);

    // Apply displacement constraints
    const pointConstraints& pcs = pointConstraints::New(pointD_.mesh());
    pcs.constrainDisplacement(pointD_, false);


    // Update displacements
    x_ = mesh().C() + D_;
    DD() = D_ - D_.oldTime();
    pointDD() = pointD_ - pointD_.oldTime();

    U_.ref() = rhoU_()/rho()();
    U_.boundaryFieldRef() =
        DD().boundaryField()/mesh().time().deltaTValue();
    rhoU_.boundaryFieldRef() = rho().boundaryField()*U_.boundaryField();

    return true;
}


tmp<vectorField> explicitTotalLagSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch mechanical property
    const scalarField impK(mechanical().impK(patchID));

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals
    const vectorField& n(mech_.n().boundaryField()[patchID]);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - n*pressure)
              - (n & (pSigma - impK*pGradD))
            )/impK
        )
    );
}


scalar explicitTotalLagSolid::CoNum() const
{
    // waveSpeed = cellWidth/deltaT
    // So, deltaT = cellWidth/waveVelocity == (1.0/deltaCoeff)/waveSpeed
    // In the current discretisation, information can move two cells per
    // time-step. This means that we use 1/(2*d) == 0.5*deltaCoeff when
    // calculating the required stable time-step
    // i.e.e deltaT = (1.0/(0.5*deltaCoeff)/waveSpeed
    // For safety, we should use a time-step smaller than this e.g. Abaqus uses
    // 1/sqrt(2)*stableTimeStep: we will default to this value

    surfaceScalarField amaxSf
    (
        fvc::interpolate
        (
            sqrt(mechanical().elasticModulus()/rho())
           /mech_.stretch()
        )*mesh().magSf()
    );

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


scalar explicitTotalLagSolid::maxCoNum() const
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
