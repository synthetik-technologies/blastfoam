/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
#include "pressureReference.H"
#include "findRefCell.H"
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
    const dictionary& dict,
    const label zoneID
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
    baseElevation_
    (
        IOobject
        (
            "hRef",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dict_.found("hRef")
      ? dimensionedScalar("hRef", dimLength, dict_)
      : dimensionedScalar("hRef", dimLength, 0.0)
    ),
    h_("h", normal_ & mesh.C()),
    zoneID_(zoneID),
    fixedPatches_(dict.lookupOrDefault("fixedPatches", wordReList()))
{
    h_ += baseElevation_ - min(h_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphereModel::~atmosphereModel()
{}

void Foam::atmosphereModel::hydrostaticInitialisation
(
    fluidBlastThermo& thermo
) const
{
    volScalarField p(thermo.p());

    const volScalarField& rho(thermo.rho());
    const fvMesh& mesh = p.mesh();

    volScalarField gh("gh", (g_ & normal_)*h_);
    surfaceScalarField ghf("ghf", fvc::interpolate(gh));

    // Set the default boundary conditions for ph_rgh
    wordList ph_rghBcs(p.boundaryField().size(), "fixedFluxPressure");
    forAll(ph_rghBcs, patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            ph_rghBcs[patchi] = "fixedValue";
        }
    }
    if (fixedPatches_.size())
    {
        Info<< "Fixing " << fixedPatches_ << endl;
        labelHashSet fixedPatchIDs(mesh.boundaryMesh().patchSet(fixedPatches_));
        forAllConstIter(labelHashSet, fixedPatchIDs, iter)
        {
            const label patchi = iter.key();
            ph_rghBcs[patchi] = "fixedValue";
        }
    }

    // Optionally read in some fields
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
        "fixedValue"
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

    volScalarField ph_rgh
    (
        IOobject
        (
            "ph_rgh",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh_,
        dimensionedScalar("ph_rgh", dimPressure, 0.0),
        ph_rghBcs
    );
    ph_rgh.rename("p");

    pressureReference pressureReference
    (
        ph_rgh,
        dict_,
        ph_rgh.needReference()
    );

    label nCorr
    (
        dict_.lookupOrDefault<label>("nHydrostaticCorrectors", 10)
    );

    // Create a simple solver dictionary
    dictionary solverDict;
    solverDict.add("solver", "PCG");
    solverDict.add("preconditioner", "DIC");
    solverDict.add("smoother", "GaussSeidel");
    solverDict.add("tolerance", 1e-6);
    solverDict.add("relTol", 0);
    solverDict.add("minIter", 1);

    for (label i=0; i<nCorr; i++)
    {
        Info<< nl << "Hydrostatic iteration " << i << endl;

        ph_rgh == p - rho*gh;

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
        ph_rghEqn.setReference
        (
            pressureReference.refCell(),
            getRefCellValue(ph_rgh, pressureReference.refCell())
        );

        ph_rghEqn.solve(solverDict);


        scalar residual = (max(ph_rgh()) - min(ph_rgh())).value();
        p = ph_rgh + rho*gh;
        if (ph_rgh.needReference())
        {
            p += dimensionedScalar
            (
                "pRef",
                p.dimensions(),
                pressureReference.refValue()
              - getRefCellValue(p, pressureReference.refCell())
            );
        }

        Info<< "Hydrostatic pressure variation "<< residual << endl;

        // Correct density and thermodynamic quantities
        p.correctBoundaryConditions();
        thermo.updateRho(p);
        thermo.he() = thermo.calce(thermo.p());
    }

    if (zoneID_ >= 0)
    {
        const cellZone& cz = mesh_.cellZones()[zoneID_];
        UIndirectList<scalar>(thermo.p(), cz) = UIndirectList<scalar>(p, cz);
        thermo.p().correctBoundaryConditions();

        // Correct density and thermodynamic quantities
        thermo.updateRho(thermo.p());
        thermo.he() = thermo.calce(thermo.p());
        thermo.correct();
    }
}

// ************************************************************************* //
