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
    groundElevation_("groundElevation", dimLength, dict_),
    h_("h", normal_ & mesh.C())
{
    h_ += groundElevation_ - min(h_);
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

    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            mesh.time().timeName(),
            mesh
        ),
        p - rho*gh - pRef,
        "fixedFluxPressure"
    );

    volVectorField U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        "noSlip"
    );

    volScalarField ph_rgh
    (
        IOobject
        (
            "ph_rgh",
            mesh
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -2, 0, 0, 0, 0), 0.0),
        "fixedFluxPressure"
    );

    label nCorr
    (
        dict_.lookupOrDefault<label>("nHydrostaticCorrectors", 5)
    );

    dictionary solverDict;
    solverDict.add("solver", "PCG");
    solverDict.add("preconditioner", "DIC");
    solverDict.add("tolerance", 1e-6);
    solverDict.add("relTol", 0);

    for (label i=0; i<nCorr; i++)
    {
        surfaceScalarField rhof("rhof", fvc::interpolate(rho));

        surfaceScalarField phig
        (
            "phig",
            -rhof*ghf*fvc::snGrad(rho)*mesh.magSf()
        );

        // Update the pressure BCs to ensure flux consistency
        constrainPressure(ph_rgh, rho, U, phig, rhof);

        fvScalarMatrix ph_rghEqn
        (
            fvm::laplacian(rhof, ph_rgh) == fvc::div(phig)
        );

        ph_rghEqn.solve(solverDict);

        p = ph_rgh + rho*gh + pRef;
        thermo.updateRho();
        thermo.calce(p);

        Info<< "Hydrostatic pressure variation "
            << (max(ph_rgh) - min(ph_rgh)).value() << endl;
    }
}

// ************************************************************************* //
