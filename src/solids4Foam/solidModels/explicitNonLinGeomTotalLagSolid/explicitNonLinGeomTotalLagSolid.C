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

#include "explicitNonLinGeomTotalLagSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "wedgeFvPatch.H"
#include "wedgeFvPatchField.H"
#include "symmetryPolyPatch.H"
#include "symmetryPlanePolyPatch.H"
#include "solidTractionFvPatchVectorField.H"
#include "fixedValueFvPatchFields.H"

#include "calculatedLinearMomentumFvPatchVectorField.H"
#include "calculatedTractionFvPatchVectorField.H"

#include "symmetricLinearMomentumFvPatchVectorField.H"
#include "symmetricTractionFvPatchVectorField.H"

#include "tractionLinearMomentumFvPatchVectorField.H"
#include "tractionTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(explicitNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, explicitNonLinGeomTotalLagSolid, dictionary
);

// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

wordList
explicitNonLinGeomTotalLagSolid::linearMomentumBoundaryTypes() const
{
    wordList bTypes
    (
        mesh().boundary().size(),
        calculatedLinearMomentumFvPatchVectorField::typeName
    );
    const volVectorField& Disp = D();
    forAll(bTypes, patchi)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchi];
        const fvPatchVectorField& pD(Disp.boundaryField()[patchi]);
        if (isA<solidTractionFvPatchVectorField>(pD))
        {
            bTypes[patchi] =
                tractionLinearMomentumFvPatchVectorField::typeName;
        }
        else if
        (
            isA<symmetryPolyPatch>(patch)
         || isA<symmetryPlanePolyPatch>(patch)
        )
        {
            bTypes[patchi] =
                symmetricLinearMomentumFvPatchVectorField::typeName;
        }
        else if (polyPatch::constraintType(patch.type()))
        {
            bTypes[patchi] = patch.type();
        }
    }

    return bTypes;
}


wordList
explicitNonLinGeomTotalLagSolid::tractionBoundaryTypes() const
{
    wordList bTypes
    (
        mesh().boundary().size(),
        calculatedTractionFvPatchVectorField::typeName
    );
    const volVectorField& Disp = D();
    forAll(bTypes, patchi)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchi];
        const fvPatchVectorField& pD(Disp.boundaryField()[patchi]);
        if (isA<solidTractionFvPatchVectorField>(pD))
        {
            bTypes[patchi] =
                tractionTractionFvPatchVectorField::typeName;
        }
        else if
        (
            isA<symmetryPolyPatch>(patch)
         || isA<symmetryPlanePolyPatch>(patch)
        )
        {
            bTypes[patchi] =
                symmetricTractionFvPatchVectorField::typeName;
        }
        else if (polyPatch::constraintType(patch.type()))
        {
            bTypes[patchi] = patch.type();
        }
    }

    return bTypes;
}


void explicitNonLinGeomTotalLagSolid::updateStress()
{
    // Update gradient of displacement
    gradD() = gradSchemes_.gradient(D());

    // Update increment of displacement
    DD() = D() - D().oldTime();

    // Update gradient of displacement increment
    gradDD() = gradD() - gradD().oldTime();

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());

    P_ = J_*(sigma() & ops_.invT(F_));

    impK_ = mechanical().impK();
    impKf_ = mechanical().impKf();
}


void explicitNonLinGeomTotalLagSolid::solveGEqns
(
    volVectorField& rhoURHS,
    volVectorField& rhoURHS1,
    const label stage
)
{
    // Compute right hand sides
    rhoURHS = fvc::div(tractionf_*mesh().magSf());

    if (angularMomentumConservation_)
    {
        volVectorField rhsRhoUAM
        (
            fvc::div((xf_ ^ tractionf_)*mesh().magSf())
        );
        am_.AMconservation(rhoURHS, rhoURHS1, rhsRhoUAM, stage);
    }

    dimensionedScalar deltaT(mesh().time().deltaT());

    // Update linear momentum
    rhoU_ += deltaT*rhoURHS;
    U() = rhoU_/rho();
    U().correctBoundaryConditions();

    // Update coordinates
    surfaceScalarField rhof(fvc::interpolate(rho()));
    Uf_ = rhoUf_/rhof;

    D() += deltaT*U();
    D().correctBoundaryConditions();
    forAll(D().boundaryField(), patchi)
    {
        const polyPatch& patch = mesh().boundaryMesh()[patchi];
        fvPatchVectorField& pD(D().boundaryFieldRef()[patchi]);
        if (!polyPatch::constraintType(patch.type()))
        {
            fvPatchVectorField& pU(U().boundaryFieldRef()[patchi]);
            const fvPatchVectorField& pDOld
            (
                D().oldTime().boundaryFieldRef()[patchi]
            );
            pU = (pD - pDOld)/deltaT.value();
            rhoU_.boundaryFieldRef()[patchi] =
                rho().boundaryField()[patchi]*U().boundaryField()[patchi];
        }
    }

    // Update deformation gradient tensor
    F_ += deltaT*fvc::div(Uf_*mesh().Sf());

    x_ = mesh().C() + D();
    xf_ = fvc::interpolate(x_);

    mech_.correctN(F_);
}


void explicitNonLinGeomTotalLagSolid::updateFluxes()
{
    J_ = det(F_);
    mech_.correct(pWaveSpeed_, sWaveSpeed_, F_);

    pWaveSpeed_ =
        sqrt(mechanical().elasticModulus()/rho())/beta_/mech_.stretch();
    sWaveSpeed_ =
        sqrt(mechanical().shearModulus()/rho())*beta_/mech_.stretch();

    updateStress();

    // Cell gradients
    const surfaceVectorField& N = mech_.N();
    const surfaceVectorField& n = mech_.n();

    volVectorField Px(ops_.decomposeTensorX(P_));
    volVectorField Py(ops_.decomposeTensorY(P_));
    volVectorField Pz(ops_.decomposeTensorZ(P_));

    volTensorField gradRhoU(gradSchemes_.gradient(rhoU_));
    volTensorField gradPx(gradSchemes_.gradient(Px));
    volTensorField gradPy(gradSchemes_.gradient(Py));
    volTensorField gradPz(gradSchemes_.gradient(Pz));

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
    tractionf_ =
        0.5*(tractionOwn + tractionNei)
      + 0.5*(stabRhoU & (rhoUNei - rhoUOwn));
    rhoUf_ =
        0.5*(rhoUOwn + rhoUNei)
      + 0.5*(stabTraction & (tractionNei - tractionOwn));

    D().correctBoundaryConditions();
    traction_.correctBoundaryConditions();
    rhoU_.correctBoundaryConditions();

    forAll(rhoU_.boundaryField(), patchi)
    {
        const fvPatch& patch = mesh().boundary()[patchi];
        const fvPatchVectorField& pRhoU(rhoU_.boundaryField()[patchi]);
        // Riemann solver for inter-processor boundaries
        if (patch.coupled())
        {
            const vectorField prhoUNei
            (
                rhoU_.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pPNei
            (
                P_.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pPxNei
            (
                Px.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pPyNei
            (
                Py.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pPzNei
            (
                Pz.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradRhoUNei
            (
                gradRhoU.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradPxNei
            (
                gradPx.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradPyNei
            (
                gradPy.boundaryField()[patchi].patchNeighbourField()
            );

            const tensorField pgradPzNei
            (
                gradPz.boundaryField()[patchi].patchNeighbourField()
            );

            const vectorField pdeltaOwn
            (
                patch.fvPatch::delta()
            );
            const vectorField pdeltaNei
            (
                patch.fvPatch::delta() - patch.delta()
            );

            const scalarField ppWaveSpeedNei
            (
                pWaveSpeed_.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField psWaveSpeedNei
            (
                sWaveSpeed_.boundaryField()[patchi].patchNeighbourField()
            );

            forAll(mesh().boundary()[patchi], facei)
            {
                const label celli = patch.faceCells()[facei];
                const vector& dOwn(pdeltaOwn[facei]);
                const vector& dNei(pdeltaNei[facei]);

                const vector rhoUOwni =
                    rhoU_[celli] + (gradRhoU[celli] & dOwn);
                const vector rhoUNeii =
                    prhoUNei[facei] + (pgradRhoUNei[facei] & dNei);

                const vector PxOwni =
                    Px[celli] + (gradPx[celli] & dOwn);
                const vector PxNeii =
                    pPxNei[facei] + (pgradPxNei[facei] & dNei);

                const vector PyOwni =
                    Py[celli] + (gradPy[celli] & dOwn);
                const vector PyNeii =
                    pPyNei[facei] + (pgradPyNei[facei] & dNei);

                const vector PzOwni =
                    Pz[celli] + (gradPz[celli] & dOwn);
                const vector PzNeii =
                    pPzNei[facei] + (pgradPzNei[facei] & dNei);

                const tensor POwni = tensor(PxOwni, PyOwni, PzOwni);
                const tensor PNeii = tensor(PxNeii, PyNeii, PzNeii);

                const scalar pWSi =
                    0.5*(pWaveSpeed_[celli] + ppWaveSpeedNei[facei]);
                const scalar sWSi =
                    0.5*(sWaveSpeed_[celli] + psWaveSpeedNei[facei]);

                const vector& Ni = N.boundaryField()[patchi][facei];
                const vector& ni = n.boundaryField()[patchi][facei];

                const tensor stabRhoUi =
                    (pWSi*ni*ni) + (sWSi*(I - (ni*ni)));
                const tensor stabTractioni =
                    ((ni*ni)/pWSi) + ((I - (ni*ni))/sWSi);

                tractionf_.boundaryFieldRef()[patchi][facei] =
                    0.5*((POwni + PNeii) & Ni)
                  + (0.5*stabRhoUi & (rhoUNeii - rhoUOwni));

                rhoUf_.boundaryFieldRef()[patchi][facei] =
                    0.5*(rhoUOwni + rhoUNeii)
                  + 0.5*(stabTractioni & ((PNeii - POwni) & Ni));
            }
        }
        else
        {
            tractionf_.boundaryFieldRef()[patchi] =
                traction_.boundaryField()[patchi];

            rhoUf_.boundaryFieldRef()[patchi] = pRhoU;
        }
    }

    // Nodal linear momentum
    volVectorField rhoUAvg
    (
        interpSchemes_.surfaceToVol(rhoUf_, pointRhoU_)
    );
    volTensorField gradRhoUAvg
    (
        gradSchemes_.localGradient(rhoUAvg, rhoUf_, pointRhoU_)
    );
    interpSchemes_.volToPoint(rhoUAvg, gradRhoUAvg, pointRhoU_);

    // Symmetric boundary patch
    pointRhoU_.correctBoundaryConditions();

    // Constrained fluxes
    rhoUf_ = interpSchemes_.pointToSurface(pointRhoU_);

    // Update coordinates
    Uf_ = rhoUf_/fvc::interpolate(rho());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitNonLinGeomTotalLagSolid::explicitNonLinGeomTotalLagSolid
(
    dynamicFvMesh& mesh
)
:
    solidModel(typeName, mesh, nonLinGeom(), incremental()),
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
    J_
    (
        IOobject
        (
            "J",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        det(F_)
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
        mesh.C()
    ),
    xf_
    (
        IOobject
        (
            "xf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.Cf()
    ),
    pointX_
    (
        IOobject
        (
            "pointX",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimLength, Zero)
    ),
    Uf_
    (
        IOobject
        (
            "Uf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U())
    ),
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
        rho()*U(),
        linearMomentumBoundaryTypes(),
        mesh.boundaryMesh().types()
    ),
    rhoUf_
    (
        IOobject
        (
            "rhoUf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(rho())*Uf_
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
    traction_
    (
        IOobject
        (
            "traction",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("0", dimensionSet(1, -1, -2, 0, 0, 0, 0), Zero),
        tractionBoundaryTypes(),
        mesh.boundaryMesh().types()
    ),
    tractionf_
    (
        IOobject
        (
            "tractionf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(traction_.dimensions(), Zero)
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
        J_*(sigma() & ops_.invT(F_))
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
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf())
{
    mechanical().setUseSolidDeformation();

    pointX_.primitiveFieldRef() = mesh.points();

    DisRequired();

    // Nodal linear momentum
    volVectorField rhoUAvg
    (
        interpSchemes_.surfaceToVol(rhoUf_, pointRhoU_)
    );
    volTensorField gradRhoUAvg
    (
        gradSchemes_.localGradient(rhoUAvg, rhoUf_, pointRhoU_)
    );
    interpSchemes_.volToPoint(rhoUAvg, gradRhoUAvg, pointRhoU_);

    // Symmetric boundary patch
    pointRhoU_.correctBoundaryConditions();

    // Constrained fluxes
    rhoUf_ = interpSchemes_.pointToSurface(pointRhoU_);

    // Update coordinates
    Uf_ = rhoUf_/fvc::interpolate(rho());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void explicitNonLinGeomTotalLagSolid::setDeltaT(Time& runTime)
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


bool explicitNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    rhoU_.oldTime();
    F_.oldTime();
    D().oldTime();

    volVectorField rhoURHS
    (
        volVectorField::New
        (
            "rhoURHS",
            mesh(),
            dimensionedVector(rhoU_.dimensions()/dimTime, Zero)
        )
    );
    volVectorField rhoURHS1(rhoURHS);

    updateFluxes();

    solveGEqns(rhoURHS, rhoURHS1, 0);

    updateFluxes();
    solveGEqns(rhoURHS, rhoURHS1, 1);

    rhoU_ = 0.5*(rhoU_.oldTime() + rhoU_);
    F_ = 0.5*(F_.oldTime() + F_);
    D() = 0.5*(D().oldTime() + D());

    // Update displacements
    DD() == D() - D().oldTime();
    mechanical().interpolate(D(), pointD());
    pointDD() = pointD() - pointD().oldTime();

        // Check energies
//         energies_.checkEnergies
//         (
//             rho(), U(), D(), DD(), sigma(), gradD(), gradDD(), pWaveSpeed_, g(),
//             0.0, impKf_
//         );

    return true;
}


tmp<vectorField> explicitNonLinGeomTotalLagSolid::tractionBoundarySnGrad
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


scalar explicitNonLinGeomTotalLagSolid::CoNum() const
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
        fvc::interpolate(pWaveSpeed_)*mesh().magSf()
    );

    // Remove wave speed from wedge boundaries
    forAll(amaxSf.boundaryField(), patchi)
    {
        if (isA<wedgeFvPatch>(mesh().boundary()[patchi]))
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


scalar explicitNonLinGeomTotalLagSolid::maxCoNum() const
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
