/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "atmosphereModel.H"
#include "fluidThermo.H"
#include "fvmLaplacian.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "fvcFlux.H"
#include "surfaceInterpolate.H"
#include "constrainPressure.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(atmosphereModel, 0);
    defineRunTimeSelectionTable(atmosphereModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphereModel::atmosphereModel
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    dict_(dict),
    mesh_(mesh),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    normal_(-g_.value()/max(mag(g_).value(), small)),
    baseElevation_("baseElevation", dimLength, dict_),
    h_("h", normal_ & mesh.C())
{
    h_ += baseElevation_ - min(h_);
}


Foam::atmosphereModel::atmosphereModel
(
    const fvMesh& mesh,
    const scalar elevation
)
:
    dict_(),
    mesh_(mesh),
    g_
    (
        IOobject
        (
            "g",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    normal_(-g_.value()/max(mag(g_).value(), small)),
    baseElevation_("baseElevation", dimLength, elevation),
    h_("h", normal_ & mesh.C())
{
    h_ += baseElevation_ - min(h_);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphereModel::~atmosphereModel()
{}

void Foam::atmosphereModel::hydrostaticInitialisation
(
    fluidBlastThermo& thermo,
    const dimensionedScalar& pRef
) const
{
    volScalarField& p(thermo.p());
    volScalarField& rho(thermo.rho());
    const fvMesh& mesh = p.mesh();

    volScalarField gh("gh", (g_ & normal_)*h_);
    surfaceScalarField ghf("ghf", fvc::interpolate(gh));

    wordList p_rghBcs(p.boundaryField().size(), "fixedFluxPressure");
    forAll(p_rghBcs, patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            p_rghBcs[patchi] = "fixedValue";
        }
    }

    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        "noSlip"
    );
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh
        ),
        fvc::flux(U)
    );

    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0),
        p_rghBcs
    );

    label nCorr
    (
        dict_.lookupOrDefault<label>("nHydrostaticCorrectors", 10)
    );
    scalar tolerance(dict_.lookupOrDefault<scalar>("tolerance", 1e-6));

    dictionary solverDict;
    solverDict.add("solver", "PCG");
    solverDict.add("preconditioner", "DIC");
    solverDict.add("tolerance", 1e-10);
    solverDict.add("relTol", 0.1);

    surfaceScalarField rhof("rhof", fvc::interpolate(rho));
    surfaceScalarField phig
    (
        "phig",
        -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
    );

    for (label i=0; i<nCorr; i++)
    {
        Info<< nl << "Hydrostatic iteration " << i << endl;

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(p_rgh, rho, U, phig, rhof);

        fvScalarMatrix ph_rghEqn
        (
            fvm::laplacian(rhof, p_rgh) == fvc::div(phig)
        );

        ph_rghEqn.solve(solverDict);

        scalar residual = (max(p_rgh) - min(p_rgh)).value();

        p = p_rgh + rho*gh + pRef;
        p.correctBoundaryConditions();

        p_rgh = p - rho*gh - pRef;

        Info<< "Hydrostatic pressure variation "<< residual << endl;

        if (residual < tolerance)
        {
            break;
        }
    }
    thermo.correct();
    p_rgh.write();
}

// ************************************************************************* //
