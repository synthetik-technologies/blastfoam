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
#include "globalInterpolatedPointPatchFields.H"
#include "globalMappedPointPatchFields.H"
#include "slipPointPatchFields.H"
#include "fixedValuePointPatchFields.H"
#include "cellMotionFvPatchFields.H"
#include "motionDiffusivity.H"

#include "fvmLaplacian.H"
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

Foam::regionSolvers::fluid::fluid(dynamicFvMesh& mesh)
:
    regionSolver(mesh),
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
    velocityFields_(1, "U")
{
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

    const dictionary& dynMeshDict = dynMesh_.dynamicMeshDict();
    if (!dynMeshDict.found("diffusivity"))
    {
        OStringStream os;
        os  << "quadratic" << token::SPACE << "inverseDistance"
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::fluid::~fluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::regionSolvers::fluid::initialiseMesh(const bool firstIter)
{
    if (firstIter)
    {
        pointDPtr_().primitiveFieldRef() = Zero;
        forAll(pointDPtr_->boundaryField(), patchi)
        {
            if
            (
                isA<valuePointPatchVectorField>
                (
                    pointDPtr_->boundaryField()[patchi]
                )
            )
            {
                dynamicCast<valuePointPatchVectorField>
                (
                    pointDPtr_->boundaryFieldRef()[patchi]
                ) == Zero;
            }
        }

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
        }
        else
        {
            pointsOldPtr_.reset(new pointField(mesh_.points()));
        }
    }

    moveMesh(true);

    if (mesh_.moving())
    {
        const_cast<surfaceScalarField&>(mesh_.phi()) == Zero;
    }

    return false;
}


void Foam::regionSolvers::fluid::initialise()
{}


bool Foam::regionSolvers::fluid::changeMesh()
{
    bool changed = regionSolver::changeMesh();
    pointsOldPtr_.reset(nullptr);
    pointsOldPtr_.set(new pointField(mesh_.points()));

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


bool Foam::regionSolvers::fluid::moveMesh(const bool finalIter)
{
    regionSolver::moveMesh(finalIter);

    pointDPtr_->oldTime();
    pointDPtr_->storePrevIter();

    // Solve point motion

    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    mesh_.movePoints(pointsOldPtr_());

    diffusivityPtr_->correct();
    pointDPtr_->boundaryFieldRef().updateCoeffs();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellDPtr_(),
            "laplacian(diffusivity,cellD)"
        )
    );


    // Update point displacement
    volPointInterpolation::New(mesh_).interpolate
    (
        cellDPtr_(),
        pointDPtr_()
    );

    if (gMax(mag(pointDPtr_->primitiveField())) < small)
    {
        return false;
    }

    if (!finalIter && mesh().relaxField(cellDPtr_->name()))
    {
        scalar f = mesh().fieldRelaxationFactor(cellDPtr_->name());
        pointDPtr_() == pointDPtr_().prevIter()*(1.0 - f) + f*pointDPtr_();
    }

    tmp<pointField> tcurPoints
    (
        pointsOldPtr_()
      + (pointDPtr_->primitiveField() - pointDPtr_->oldTime().primitiveField())
    );

    twoDPointCorrector::New(mesh_).correctPoints(tcurPoints.ref());
    mesh_.movePoints(tcurPoints);

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
    return true;
}

// ************************************************************************* //
