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
    gradDD() = gradD() - gradD().oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
    impK_ = mechanical().impK();

    P_ = mech_.J()*(sigma() & ops_.invT(F_));
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
        am_.AMconservation(x_, rhoU_, rhoURHS, rhsRhoUAM, stage);
    }
//     vector validD((mesh().geometricD() + labelVector::one)/2);
//     rhoURHS = cmptMultiply(validD, rhoURHS);

    // Update coordinates
    volVectorField prevD(D_);
    U_ = rhoU_/rho();
    D_ += deltaT*U_;
    D_.correctBoundaryConditions();

    x_ = mesh().C() + D_;

    // Update linear momentum
    rhoU_ += deltaT*rhoURHS;
    rhoU_.correctBoundaryConditions();

    pointVectorField prevPointD(pointD_);
    mechanical().interpolate(D_, pointD_);
    Df_ = interpSchemes_.pointToSurface(pointD_);

    // Update increment of displacement
    DD() = D_ - DD();
    pointDD() = pointD_ - pointDD();

    // Update deformation gradient tensor
    F_ += deltaT*fvc::div(rhoUC_/fvc::interpolate(rho_)*mesh().Sf());
}


void explicitTotalLagSolid::updateFluxes()
{
    mech_.correct(pWaveSpeed_, sWaveSpeed_);

    pWaveSpeed_ =
        sqrt(mechanical().elasticModulus()/rho_)/beta_/mech_.stretch();
    sWaveSpeed_ =
        sqrt(mechanical().shearModulus()/rho_)*beta_/mech_.stretch();

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
        const vectorField pn(n.boundaryField()[patchi]);
        const vectorField pN(N.boundaryField()[patchi]);

        // Riemann solver for inter-processor boundaries
        if (pD.fixesValue() || U_.boundaryField()[patchi].fixesValue())
        {
            prhoUC[patchi] =
                rho_.boundaryField()[patchi]
               *DD().boundaryField()[patchi]
              /mesh().time().deltaT0().value();

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
            prhoUC[patchi] =
                rho_.boundaryField()[patchi]
               *(
                   DD().boundaryField()[patchi]
                )/mesh().time().deltaT().value();

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
            prhoUC[patchi] = rhoU_.boundaryField()[patchi];
            ptractionC[patchi] =  P_.boundaryField()[patchi] & pn;
        }
    }

    // Average linear momentum
    volVectorField rhoUAvg
    (
        interpSchemes_.surfaceToVol(rhoUC_, pointRhoU_)
    );
    volTensorField gradRhoUAvg
    (
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
//         else if
//         (
//             D_.boundaryField()[patchi].fixesValue()
//          || U_.boundaryField()[patchi].fixesValue()
//         )
//         {
//             dynamicCast<valuePointPatchVectorField&>(ppointRhoU[patchi]) =
//                 globalPatches()[p].faceToPoint(prhoUC[patchi]);
//         }
//         else if (!polyPatch::constraintType(p.type()))
//         {
//             ppointRhoU[patchi].setInInternalField
//             (
//                 pointRhoU_.primitiveFieldRef(),
//                 globalPatches()[p].faceToPoint(prhoUC[patchi])()
//             );
//         }
    }
    pointRhoU_.correctBoundaryConditions();
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
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
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
        sqrt(mechanical().elasticModulus()/rho_)/beta_/mech_.stretch()
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
        sqrt(mechanical().shearModulus()/rho_)*beta_/mech_.stretch()
    ),
    rhoUOwn_("rhoUOwn", rhoUC_),
    rhoUNei_("rhoUNei", rhoUC_),
    tractionOwn_("rhoUOwn", tractionC_),
    tractionNei_("rhoUNei", tractionC_),
    JSTScaleFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("JSTScaleFactor", 0.0)
    ),
    energies_(mesh, solidModelDict()),
    impK_(mechanical().impK())
{
    // Tell the mechanical laws to lookup deformation fields
    mechanical().setUseSolidDeformation();

    DisRequired();
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

    D_.oldTime();
    F_.oldTime();
    rhoU_.oldTime();
    pointD_.oldTime();

    do
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

        updateFluxes();
        solveGEqns(rhoURHS, 0);

        updateFluxes();
        solveGEqns(rhoURHS, 1);

        D_ = 0.5*(D_.oldTime() + D_);
        D_.correctBoundaryConditions();
        x_ = mesh().C() + D_;

        F_ = 0.5*(F_.oldTime() + F_);
        F_.correctBoundaryConditions();

        rhoU_ = 0.5*(rhoU_.oldTime() + rhoU_);
        rhoU_.correctBoundaryConditions();

        mechanical().interpolate(D_, pointD_);
        Df_ = fvc::interpolate(D_);

        // Update increment of displacements
        DD() = D_ - D_.oldTime();
        pointDD() = pointD_ - pointD_.oldTime();
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
