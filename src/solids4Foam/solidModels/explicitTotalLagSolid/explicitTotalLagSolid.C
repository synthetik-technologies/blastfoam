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
    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
    impK_ = mechanical().impK();

    P_ = mech_.J()*(sigma() & mech_.invF().T());
}


void explicitTotalLagSolid::solveGEqns
(
    volVectorField& rhoURHS,
    const label stage
)
{
    dimensionedScalar deltaT(mesh().time().deltaT());

    // Compute right hand sides
    rhoURHS = fvc::surfaceIntegrate(tractionC_*mesh().magSf());
    if (JSTScaleFactor_ > 0)
    {
        rhoURHS -=
            JSTScaleFactor_*fvc::laplacian
            (
                mesh().magSf(),
                fvc::laplacian
                (
                    deltaT*mechanical().impKf(),
                    U_,
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
        am_.AMconservation(x_, rhoU_, rhoURHS, rhsRhoUAM, stage);
    }

    surfaceScalarField rhof(fvc::interpolate(rho_));

    // Update coordinates
    {
        volVectorField prevD(D_);
        volTensorField prevGradD(gradD_);
        pointVectorField prevPointD(pointD_);

        D_.ref() += U_()*deltaT;
        D_.correctBoundaryConditions();

        Df_ += deltaT*rhoUC_/rhof;
        pointD_ += deltaT*pointRhoU_/mechanical().volToPoint().interpolate(rho_);;

//         mechanical().interpolate(D_, pointD_);
        mechanical().grad(D_, gradD_);

        DD_ = D_ - prevD;
        gradDD_ = gradD_ - prevGradD;
        pointDD_ = pointD_ - prevPointD;
    }

    // Face displacement
    Df_ = fvc::interpolate(D_);

    // Material positions
    x_ = mesh().C() + D_;

    // Update linear momentum
    rhoU_.ref() += deltaT*rhoURHS();
    U_.ref() = rhoU_()/rho_();

    U_.boundaryFieldRef() = DD_.boundaryField()/deltaT.value();
    rhoU_.boundaryFieldRef() = U_.boundaryField()*rho_.boundaryField();

    // Update deformation gradient tensor
    F_ +=
        deltaT
       *fvc::surfaceIntegrate(rhoUC_/rhof*mesh().Sf());
//     F_ = I + gradD_.T();

    // Update deformation quantities
    mech_.correctDeformation();
}


void explicitTotalLagSolid::updateFluxes()
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

    volVectorField Px(ops_.decomposeTensorX(P_));
    volVectorField Py(ops_.decomposeTensorY(P_));
    volVectorField Pz(ops_.decomposeTensorZ(P_));

    volTensorField gradRhoU(fvc::grad(rhoU_));
    volTensorField gradPx(fvc::grad(Px));
    volTensorField gradPy(fvc::grad(Py));
    volTensorField gradPz(fvc::grad(Pz));

    // Reconstruction
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


    gradSchemes_.reconstruct(rhoU_, gradRhoU, rhoUOwn_, rhoUNei_);
    gradSchemes_.reconstruct(P_, gradPx, gradPy, gradPz, POwn, PNei);
    tractionOwn_ = POwn & N;
    tractionNei_ = PNei & N;

    const surfaceTensorField& stabRhoU(mech_.stabRhoU());
    const surfaceTensorField& stabTraction(mech_.stabTraction());

    // Acoustic Riemann solver
    tractionC_ =
        0.5*(tractionOwn_ + tractionNei_)
      + 0.5*(stabRhoU & (rhoUNei_ - rhoUOwn_));
    rhoUC_ =
        0.5*(rhoUOwn_ + rhoUNei_)
      + 0.5*(stabTraction & (tractionNei_ - tractionOwn_));

    surfaceVectorField::Boundary& prhoUC(rhoUC_.boundaryFieldRef());
    surfaceVectorField::Boundary& ptractionC(tractionC_.boundaryFieldRef());
    forAll(ptractionC, patchi)
    {
        const polyPatch& p = mesh().boundaryMesh()[patchi];
        const fvPatch& patch = mesh().boundary()[patchi];
        const fvPatchField<vector>& pD(D_.boundaryField()[patchi]);
        const fvPatchField<vector>& pRhoU(rhoU_.boundaryField()[patchi]);
        const vectorField pn(n.boundaryField()[patchi]);
        const vectorField pN(N.boundaryField()[patchi]);

        // Riemann solver for inter-processor boundaries
        if (pD.fixesValue() || U_.boundaryField()[patchi].fixesValue())
        {
            prhoUC[patchi] = pRhoU[patchi];

            ptractionC[patchi] =
                tractionOwn_.boundaryField()[patchi]
              + (
                    stabRhoU.boundaryField()[patchi]
                  & (
                        prhoUC[patchi]
                      - rhoUOwn_.boundaryField()[patchi]
                    )
                );
        }
        else if (isA<solidTractionFvPatchVectorField>(pD))
        {
            const solidTractionFvPatchVectorField& stD =
                dynamicCast<const solidTractionFvPatchVectorField>(pD);
            vectorField tp(stD.traction() - stD.pressure()*pn);

            prhoUC[patchi] =
                rhoUOwn_.boundaryField()[patchi]
              + (
                    stabTraction.boundaryField()[patchi]
                  & (tp - tractionOwn_.boundaryField()[patchi])
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
                    rhoUOwn_.boundaryField()[patchi]
                  - tractionOwn_.boundaryField()[patchi]
                   /sWaveSpeed_.boundaryField()[patchi]
                );
            ptractionC[patchi] =
                (pn*pn)
              & (
                    tractionOwn_.boundaryField()[patchi]
                  - pWaveSpeed_.boundaryField()[patchi]
                   *rhoUOwn_.boundaryField()[patchi]
                );
        }
        else if (!polyPatch::constraintType(p.type()))
        {
            prhoUC[patchi] = pRhoU[patchi];

            ptractionC[patchi] =
                tractionOwn_.boundaryField()[patchi]
              + (
                    stabRhoU.boundaryField()[patchi]
                  & (
                        prhoUC[patchi]
                      - rhoUOwn_.boundaryField()[patchi]
                    )
                );
        }
        else
        {
            prhoUC[patchi] = pRhoU[patchi];
            ptractionC[patchi] =  P_.boundaryField()[patchi] & pn;
        }
    }

    // Average linear momentum
    volVectorField rhoUAvg
    (
//         fvc::average(rhoUC_)
        interpSchemes_.surfaceToVol(rhoUC_, pointRhoU_)
    );
    volTensorField gradRhoUAvg
    (
//         fvc::grad(rhoUAvg)
        gradSchemes_.localGradient(rhoUAvg, rhoUC_, pointRhoU_)
    );

    interpSchemes_.volToPoint(rhoUAvg, gradRhoUAvg, pointRhoU_, true);

    forAll(ptractionC, patchi)
    {
        const polyPatch& p = mesh().boundaryMesh()[patchi];
        const vectorField pN(N.boundaryField()[patchi]);
        if
        (
            isA<symmetryPolyPatch>(p)
         || isA<symmetryPlanePolyPatch>(p)
        )
        {
            forAll(p, fi)
            {
                const label& facei = p.start() + fi;
                forAll(mesh().faces()[facei], ni)
                {
                    const label& nodei = mesh().faces()[facei][ni];
                    pointRhoU_[nodei] =
                        (
                            tensor::I - (pN[facei]*pN[facei])
                        ) & pointRhoU_[nodei];
                }
            }
        }
    }
    pointRhoU_.correctBoundaryConditions();

    rhoUC_.ref() = interpSchemes_.pointToSurface(pointRhoU_)()();
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
        dimensionedVector("0", rhoU_.dimensions(), Zero),
        pointDBoundaryTypes(D_)
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
    JSTScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("JSTScaleFactor", 0.0)
    ),
    energies_(mesh, solidModelDict()),
    impK_(mechanical().impK())
{
    DisRequired();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool explicitTotalLagSolid::evolve()
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
        rhoU_.ref() = 0.5*(rhoU_.oldTime()() + rhoU_());
        U_.ref() = rhoU_()/rho_();

        D_ = 0.5*(D_.oldTime() + D_);
        Df_ = 0.5*(Df_.oldTime() + Df_);
        pointD_ = 0.5*(pointD_.oldTime() + pointD_);

        D_.correctBoundaryConditions();
        DD_ = D_ - D_.oldTime();
        pointDD_ = pointD_ - pointD_.oldTime();

        // Update coordinates
        x_ = mesh().C() + D_;

//         D_.ref() =  D_.oldTime() + deltaT*U_();
//         D_.correctBoundaryConditions();
//         mechanical().interpolate(D_, pointD_);
//
//         // Displacements
//         DD_ = D_ - D_.oldTime();
//         Df_ = fvc::interpolate(D_);
//         pointDD_ = pointD_ - pointD_.oldTime();
//
//         // Material positions
//         x_ = mesh().C() + D_;
//
        U_.boundaryFieldRef() = DD_.boundaryField()/deltaT.value();
        rhoU_.boundaryFieldRef() = U_.boundaryField()*rho_.boundaryField();

        // Update gradient of displacement increment
        mechanical().grad(D_, gradD_);

        // Update the gradient of total displacement
        gradDD_ = gradD_ - gradD_.oldTime();
//
//         // Update deformation gradient tensor
//         F_ = 0.5*(F_.oldTime() + F_);//I + gradD_.T();

        // Update deformation quantities
        mech_.correctDeformation(true);
    } while (mesh().update());

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
    const scalarField& impK(impK_.boundaryField()[patch.index()]);

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
