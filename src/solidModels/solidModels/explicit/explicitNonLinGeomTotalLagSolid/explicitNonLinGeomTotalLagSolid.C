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
#include "extrapolatedCalculatedFvPatchFields.H"
#include "wedgeFvPatch.H"

#include "operations.H"
#include "mechanics.H"
#include "gradientSchemes.H"


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

void explicitNonLinGeomTotalLagSolid::updateStress()
{
    totalLagSolid<totalDispSolid>::update(true);

    waveSpeed_ = sqrt(this->impKf_/fvc::interpolate(rho()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

explicitNonLinGeomTotalLagSolid::explicitNonLinGeomTotalLagSolid
(
    dynamicFvMesh& mesh
)
:
    totalLagSolid<totalDispSolid>(typeName, mesh),
    waveSpeed_
    (
        IOobject
        (
            "waveSpeed",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(Foam::sqrt(this->impK_/rho()))
    ),
    energies_(mesh, solidModelDict()),
    a_
    (
        IOobject
        (
            "a",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity/dimTime, Zero),
        extrapolatedCalculatedFvPatchScalarField::typeName
    ),
    relaxation_(solidModelDict().optionalSubDict("relaxation"))
{
    D().correctBoundaryConditions();
    // Update stress
    updateStress();

    // Update initial acceleration
    a_ = fvc::div(sigma(), "div(sigma)")/rho();
    a_.correctBoundaryConditions();
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
        1.0/max(mesh().surfaceInterpolation::deltaCoeffs()*waveSpeed_).value();

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

    relaxation_.read(solidModelDict().optionalSubDict("relaxation"));

    enforceLinear() = false;

    // Mesh update loop
    do
    {
        Info<< "Solving the momentum equation for D" << endl;

        // Central difference scheme
        const dimensionedScalar& deltaT = time().deltaT();
        const dimensionedScalar deltaT01(0.5*(deltaT + time().deltaT0()));

        // Compute the velocity
        // Note: this is the velocity at the middle of the time-step
        U() = U().oldTime() + deltaT01*a_.oldTime();

        // Compute displacement
        D() = D().oldTime() + deltaT*U();

        // Enforce any cell displacements
        if (setCellDisps().cellIDs().size())
        {
            UIndirectList<vector>
            (
                D(),
                setCellDisps().cellIDs()
            ) = setCellDisps().cellDisps();
        }

        // Enforce boundary conditions on the displacement field
        D().correctBoundaryConditions();

        U() = (D() - D().oldTime())/deltaT;

        // Update the stress field based on the latest D field
        updateStress();

        // Compute acceleration
        // Note the inclusion of a linear bulk viscosity pressure term to
        // dissipate high frequency energies, and a Rhie-Chow term to
        // avoid checker-boarding
        tmp<volVectorField> stab
        (
            stabilisation().stabilisation
            (
                U(),
                fvc::grad(U())(),
               (deltaT01*impKf_)()
            )
        );

        a_ =
            (
                fvc::div(this->P(), "div(sigma)")
              + fvc::div
                (
                    mesh().Sf()*energies_.viscousPressure
                    (
                        rho(), waveSpeed_, gradD()
                    )
                )

                // This corresponds to Lax–Friedrichs smoothing
              + stab()
            )/rho()
          + g();
        a_.correctBoundaryConditions();

        // Check energies
        energies_.checkEnergies
        (
            rho(),
            U(),
            D(),
            DD(),
            this->sigma(),
            gradD(),
            gradDD(),
            stab(),
            g()
        );

    } while (mesh().update());

    relaxation_.relax(U(), rho());

    return true;
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

    surfaceScalarField amaxSf(waveSpeed_*mesh().magSf());

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
