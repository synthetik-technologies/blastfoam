/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "fluidRegionSolver.H"
#include "surfaceFields.H"

#include "wallPolyPatch.H"
#include "calculatedPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "globalInterpolatedPointPatchFields.H"
#include "globalMappedPointPatchFields.H"
#include "slipPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "cellMotionFvPatchFields.H"
#include "motionDiffusivity.H"

#include "fvm.H"
#include "fvc.H"
#include "twoDPointCorrector.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(fluid, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::fluid::fluid
(
    dynamicFvMesh& mesh,
    const regionSolverList& regions
)
:
    regionSolver(mesh, regions),
    devRhoReff_
    (
        IOobject
        (
            "devRhoReff",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor
        (
            "Zero",
            dimensionSet(1,1,-2,0,0,0,0),
            Zero
        )
    ),
    relaxation_(displacementRelaxation::New(mesh_, regions_.solutionControls())),
    velocityFields_(1, "U"),
    tolerance_(-great),
    relTol_(-great)
{
    this->readControls("pointD", tolerance_, relTol_);

    const pointMesh& pMesh = pointMesh::New(mesh_);
    pointDPtr_.set
    (
        new pointVectorField
        (
            IOobject
            (
                "pointD",
                mesh_.time().timeName(),
                mesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector(dimLength, Zero),
            fixedValuePointPatchVectorField::typeName
        )
    );

    if (!pointDPtr_->headerOk())
    {
        pointVectorField::Boundary& bpointD =
            pointDPtr_->boundaryFieldRef();
        forAll(mesh_.boundary(), patchi)
        {
            const pointPatch& p = pMesh.boundary()[patchi];
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            if (globalBoundary_.isCoupled(pp))
            {
                bpointD.set
                (
                    patchi,
                    new globalMappedPointPatchVectorField
                    (
                        p,
                        pointDPtr_(),
                        "pointD"
                    )
                );
            }
            else if (p.type() == fvPatch::typeName)
            {
                bpointD.set
                (
                    patchi,
                    new slipPointPatchVectorField
                    (
                        p,
                        pointDPtr_()
                    )
                );
            }
        }
    }
    pointDPtr_->storeOldTimes();

    const pointVectorField::Boundary& bpointD =
            pointDPtr_->boundaryField();
    wordList cellDBCs(bpointD.types());
    forAll(cellDBCs, patchi)
    {
        if (isA<fixedValuePointPatchVectorField>(bpointD[patchi]))
        {
            cellDBCs[patchi] =
                cellMotionFvPatchVectorField::typeName;
        }

        if (debug)
        {
            Pout<< "Patch:" << mesh_.boundary()[patchi].patch().name()
                << " pointType:" << bpointD.types()[patchi]
                << " cellType:" << cellDBCs[patchi] << endl;
        }
    }
    cellDPtr_.set
    (
        new volVectorField
        (
            IOobject
            (
                "cellD",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh_,
            dimensionedVector(dimLength, Zero),
            cellDBCs
        )
    );
    cellDPtr_->oldTime();

    const dictionary& dynMeshDict = dynMesh_.dynamicMeshDict();
    if (!dynMeshDict.found("diffusivity"))
    {
        OStringStream os;
        os  <<  token::SPACE << "inverseDistance"
            << token::SPACE << token::BEGIN_LIST;

        forAll(mesh_.boundary(), patchi)
        {
            const polyPatch& pp = mesh_.boundaryMesh()[patchi];
            if (globalBoundary_.isCoupled(pp))
            {
                os << pp.name() << token::SPACE;
            }
        }
        os  << token::END_LIST;
        const_cast<dictionary&>(dynMeshDict).set("diffusivity", word(os.str()));
    }
    diffusivityPtr_ =
        motionDiffusivity::New(mesh_, dynMeshDict.lookup("diffusivity"));

    pointsOldPtr_.reset(new pointField(mesh_.points()));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::fluid::~fluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::fluid::initialiseMesh(const IterType iter)
{
    if (iter == FIRST_ITER)
    {
        if (mesh_.pointsInstance() != mesh_.facesInstance())
        {
            pointsOldPtr_.reset
            (
                new pointField
                (
                    pointIOField
                    (
                        IOobject
                        (
                            "points",
                            mesh_.facesInstance(),
                            polyMesh::meshSubDir,
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            false
                        )
                    )
                )
            );
            const_cast<pointField&>(mesh_.oldPoints()) = pointsOldPtr_();
            mesh_.movePoints(pointsOldPtr_());


            //- If the mesh is stored in constant, no motion has occurred
            //  So check if there is a version of pointD there
            IOobject pointDIO
            (
                pointDPtr_->name(),
                mesh_.facesInstance(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            );
            if (pointDIO.typeHeaderOk<pointVectorField>(true))
            {
                pointDPtr_() ==
                    pointVectorField
                    (
                        pointDIO,
                        pointMesh::New(mesh_)
                    );
                pointDPtr_->instance() = mesh_.time().timeName();
                pointDPtr_->storeOldTimes();
            }
            pointDPtr_() == Zero;
        }
    }

    moveMesh(FINAL_ITER);

    relaxation_->clear();

    if (mesh_.moving())
    {
        const_cast<surfaceScalarField&>(mesh_.phi()) == Zero;
        forAll(velocityFields_, i)
        {
            if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
            {
                mesh_.lookupObjectRef<volVectorField>
                (
                    velocityFields_[i]
                ).correctBoundaryConditions();
            }
        }
    }

    if (iter == FINAL_ITER)
    {
        pointDPtr_->write();
        mesh_.lookupObject<pointIOField>("points").write();
    }
}


void Foam::regionSolvers::fluid::initialise()
{
    // moveMesh(FINAL_ITER);
    // moveMesh(FINAL_ITER);
}


bool Foam::regionSolvers::fluid::changeMesh()
{
    bool changed = regionSolver::changeMesh();
    pointsOldPtr_.reset(new pointField(mesh_.points()));

    if (changed)
    {
        diffusivityPtr_.reset(nullptr);
        diffusivityPtr_ =
            motionDiffusivity::New
            (
                mesh_,
                dynMesh_.dynamicMeshDict().lookup("diffusivity")
            );
    }
    return changed;
}


bool Foam::regionSolvers::fluid::moveMesh(const IterType iter)
{
    // return true;
    regionSolver::moveMesh(iter);

    pointDPtr_->oldTime();
    pointDPtr_->storePrevIter();

    // Solve point motion

    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    mesh_.movePoints(pointsOldPtr_());

    diffusivityPtr_->correct();
    pointDPtr_->correctBoundaryConditions();

    if (iter != FINAL_ITER)
    {
        relaxation_->relax(regions_.iterNo(), pointDPtr_());
    }
    else
    {
        relaxation_->updateError(pointDPtr_());
    }

    //- Save boundary values
    pointVectorField::Boundary bpRelaxed
    (
        pointDPtr_->mesh().boundary(),
        pointDPtr_(),
        valuePointPatchVectorField::typeName
    );
    bpRelaxed = pointDPtr_->boundaryField();
    cellDPtr_->correctBoundaryConditions();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellDPtr_(),
            "laplacian(diffusivity,cellD)"
        )
    );

    // surfaceScalarField Df(diffusivityPtr_->operator()());
    // volTensorField gradCd("gradCd", fvc::grad(cellDPtr_()));
    //
    // Foam::solve
    // (
    //     fvm::laplacian
    //     (
    //         2*Df,
    //         cellDPtr_(),
    //         "laplacian(diffusivity,cellD)"
    //     )
    //
    //   + fvc::div
    //     (
    //         Df
    //        *(
    //            fvc::dotInterpolate
    //            (
    //                mesh_.Sf(),
    //                gradCd.T() - gradCd
    //            )
    //
    //            // Solid-body rotation "lambda" term
    //          - mesh_.Sf()*fvc::interpolate(tr(gradCd))
    //         )
    //     )
    //
    //
    //   // - fvc::laplacian
    //   //   (
    //   //       2*Df,
    //   //       cellDPtr_(),
    //   //       "laplacian(diffusivity,cellD)"
    //   //   )
    //   // + fvc::div
    //   //   (
    //   //       Df
    //   //      *(
    //   //          fvc::dotInterpolate
    //   //          (
    //   //              mesh_.Sf(),
    //   //              gradCd + gradCd.T()
    //   //          )
    //   //          // Solid-body rotation "lambda" term
    //   //        - mesh_.Sf()*fvc::interpolate(tr(gradCd))
    //   //      )
    //   //   )
    //
    // );


    // Update point displacement
    {
        volPointInterpolation::New(mesh_).interpolate
        (
            cellDPtr_(),
            pointDPtr_()
        );

        // Copy non-fixed values
        pointVectorField::Boundary& bp = pointDPtr_->boundaryFieldRef();
        forAll(bp, patchi)
        {
            if (isA<valuePointPatchVectorField>(bp[patchi]))
            {
                valuePointPatchVectorField& pp =
                    dynamicCast<valuePointPatchVectorField>(bp[patchi]);
                const valuePointPatchVectorField& ppRelaxed =
                    dynamicCast<const valuePointPatchVectorField>
                    (
                        bpRelaxed[patchi]
                    );
                pp == ppRelaxed;
                pp.setInInternalField(pointDPtr_(), pp);
            }
        }
    }

    tmp<pointField> tcurPoints
    (
        pointsOldPtr_()
      + (pointDPtr_->primitiveField() - pointDPtr_->oldTime().primitiveField())
    );
    twoDPointCorrector::New(mesh_).correctPoints(tcurPoints.ref());

    pointConstraints::syncUntransformedData
    (
        mesh_,
        tcurPoints.ref(),
        plusEqOp<vector>()
    );
    {
        scalarField one(tcurPoints().size(), 1.0);
        pointConstraints::syncUntransformedData
        (
            mesh_,
            one,
            plusEqOp<scalar>()
        );
        tcurPoints.ref() /= one;
    }
    mesh_.movePoints(tcurPoints());

    if (debug)
    {
        Info<<"Max mesh boundary velocity = "
            << gMax(mag(mesh_.phi().boundaryField()/mesh_.magSf().boundaryField())) <<endl;
        forAll(mesh_.boundary(), patchi)
        {
            if (!mesh_.boundary()[patchi].coupled())
            {
                Info<<mesh_.boundary()[patchi].name()<<": "
                    <<gMax(mag(mesh_.phi().boundaryField()[patchi]/mesh_.magSf().boundaryField()[patchi]))<<endl;
            }
        }
    }
    Info<<"Displacement error (abs/rel) = "
        << relaxation_->error() << ", "
        << relaxation_->relError() <<endl;

    forAll(velocityFields_, i)
    {
        if (mesh_.foundObject<volVectorField>(velocityFields_[i]))
        {
            mesh_.lookupObjectRef<volVectorField>
            (
                velocityFields_[i]
            ).correctBoundaryConditions();
        }
    }

    return gMax(mag(pointDPtr_->primitiveField())) > small;
}


void Foam::regionSolvers::fluid::clear()
{
    relaxation_->clear();
}


// ************************************************************************* //
