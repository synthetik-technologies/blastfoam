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
#include "noSlipFvPatchVectorField.H"

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
    const word& zoneName
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
    hRef_
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
    zoneName_(zoneName)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphereModel::~atmosphereModel()
{}

void Foam::atmosphereModel::hydrostaticInitialisation
(
    fluidBlastThermo& thermo
) const
{
    volScalarField p(thermo.p());

    const volScalarField& rho = thermo.rho();
    const fvMesh& mesh = p.mesh();

    volScalarField gh("gh", (g_ & mesh_.C()) + mag(g_)*hRef_);
    surfaceScalarField ghf("ghf", (g_ & mesh_.Cf()) + mag(g_)*hRef_);

    // Set the default boundary conditions for ph_rgh
    wordList ph_rghBcs(p.boundaryField().size(), "fixedFluxPressure");
    forAll(ph_rghBcs, patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            ph_rghBcs[patchi] = "fixedValue";
        }
    }
    if (dict_.found("fixedPatches"))
    {
        wordReList fixedPatches(dict_.lookup("fixedPatches"));

        Info<< "Fixing " << fixedPatches << endl;
        labelHashSet fixedPatchIDs(mesh.boundaryMesh().patchSet(fixedPatches));
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
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        noSlipFvPatchVectorField::typeName
    );
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            mesh.time().name(),
            mesh
        ),
        fvc::flux(U)
    );

    volScalarField ph_rgh
    (
        IOobject
        (
            "ph_rgh",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        p - rho*gh,
        ph_rghBcs
    );
    volScalarField T0(thermo.T());

    pressureReference pressureReference
    (
        p,
        dict_,
        ph_rgh.needReference()
    );

    label nCorr
    (
        dict_.lookupOrDefault<label>("nHydrostaticCorrectors", 10)
    );
    bool correctRho
    (
        dict_.lookupOrDefault<bool>("correctRho", true)
    );
    scalar tolerance(dict_.lookupOrDefault<scalar>("tolerance", 1e-6));
    scalar relTol(dict_.lookupOrDefault<scalar>("relTol", 1e-6));

    // Create a simple solver dictionary
    dictionary solverDict;
    solverDict.add("solver", "PCG");
    solverDict.add("preconditioner", "DIC");
    solverDict.add("smoother", "GaussSeidel");
    solverDict.add("tolerance", 1e-8);
    solverDict.add("relTol", 0);
    solverDict.add("minIter", 1);

    scalar residualOld = great;
    scalar error = great;
    label iter = 0;
    for (iter = 0; iter < nCorr; iter++)
    {
        Info<< nl << "Hydrostatic iteration " << iter << endl;

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

        // Correct density and thermodynamic quantities
        p.correctBoundaryConditions();
        if (correctRho)
        {
            thermo.T() = T0;
            thermo.updateRho(p);
        }
        thermo.he() = thermo.calce(p);
        thermo.update();

        Info<< "Hydrostatic pressure variation "<< residual << endl;
        if (iter > 0)
        {
            Info<< "Change in hydrostatic variation "
                << residual - residualOld << endl;

            error = mag(residual - residualOld);
            if (error < tolerance || error/residualOld < relTol)
            {
                break;
            }
        }
        residualOld = residual;
    }

    if (error < tolerance || error/residualOld < relTol)
    {
        Info<< nl
            << "Converged hydrostatic pressure in " << iter
            << " iterations" << nl << endl;
    }
    else
    {
        Info<< nl << "Did not converge hydrostatic pressure" << nl << endl;
    }

    if (!zoneName_.empty())
    {
        const cellZone& cz = mesh_.cellZones()[zoneName_];
        UIndirectList<scalar>(thermo.p(), cz) = UIndirectList<scalar>(p, cz);
        thermo.p().correctBoundaryConditions();

        // Correct density and thermodynamic quantities
        if (correctRho)
        {
            thermo.T() = T0;
        }
        thermo.updateRho(thermo.p());
        thermo.he() = thermo.calce(thermo.p());
        thermo.correct();
    }
    else
    {
        thermo.p() = p;
        thermo.p().correctBoundaryConditions();
    }
}

// ************************************************************************* //
