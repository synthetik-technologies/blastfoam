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
#include "syncTools.H"

#include "fvm.H"
#include "fvc.H"
#include "twoDPointCorrector.H"
#include "volPointInterpolation.H"
#include "locationMapper.H"


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
    fvMesh& mesh,
    const regionSolverList& regions
)
:
    regionSolver(mesh, regions),
    devRhoReff_
    (
        IOobject
        (
            "devRhoReff",
            runTime_.name(),
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
    velocityFields_(1, "U"),
    moving_(true)
{
    if (moving_)
    {
        const pointMesh& pMesh = pointMesh::New(mesh_);
        pointDPtr_.set
        (
            new pointVectorField
            (
                IOobject
                (
                    "pointD",
                    mesh_.time().name(),
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
                            "pointD",
                            false
                        )
                    );
                }
                else if (isA<wallPolyPatch>(pp))
                {
                    bpointD.set
                    (
                        patchi,
//                         new slipPointPatchVectorField
                        new fixedValuePointPatchVectorField
                        (
                            p,
                            pointDPtr_()
                        )
                    );
                    bpointD[patchi] == Zero;
                }
            }
        }
        pointDPtr_->oldTime();

        const pointVectorField::Boundary& bpointD =
                pointDPtr_->boundaryField();
        wordList cellDBCs(bpointD.types());
        forAll(cellDBCs, patchi)
        {
            if (isA<valuePointPatchVectorField>(bpointD[patchi]))
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
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector(dimLength, Zero),
                cellDBCs
            )
        );

        IOdictionary dynMeshDict
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        // if
        // (
        //     dynMeshDict.isDict("mover")
        //  && !dynMeshDict.subDict("mover").found("diffusivity")
        // )
        // {
        //     OStringStream os;
        //     os  <<  token::SPACE << "inverseDistance"
        //         << token::SPACE << token::BEGIN_LIST;
        //
        //     forAll(mesh_.boundary(), patchi)
        //     {
        //         const polyPatch& pp = mesh_.boundaryMesh()[patchi];
        //         if (globalBoundary_.isCoupled(pp))
        //         {
        //             os << pp.name() << token::SPACE;
        //         }
        //     }
        //     os  << token::END_LIST;
        //     const_cast<dictionary&>(dynMeshDict).set
        //     (
        //         "diffusivity",
        //         word(os.str())
        //     );
        // }
        diffusivityPtr_ =
            motionDiffusivity::New
            (
                mesh_,
                dynMeshDict.subDict("mover").lookup("diffusivity")
            );

        // Add point displacement as a relaxation field
        accelerationSchemes_.addField(pointDPtr_());

        {
            typeIOobject<pointIOField> points0IO
            (
                "points0",
                mesh_.facesInstance(),
                polyMesh::meshSubDir,
                mesh_,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE,
                false
            );

            if (points0IO.headerOk())
            {
                // Points0 written to a time folder
                points0Ptr_.set(new pointIOField(points0IO));
            }
            else
            {
                points0IO.rename("points");

                // Return copy of original mesh points
                points0Ptr_.set(new pointIOField(points0IO));
                points0Ptr_->rename("points0");
            }

            points0Ptr_->checkIn();

            locationMapper::NewRef(mesh_).addInterpolatedField
            (
                points0Ptr_->name()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::fluid::~fluid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::fluid::initialiseMesh(const IterType iter)
{
    if (moving_)
    {
        if (iter == FIRST_ITER)
        {
            if (mesh_.pointsInstance() != mesh_.facesInstance())
            {
                //- Attempt to read points0 from the lastest written mesh
                typeIOobject<pointIOField> points0IO
                (
                    "points0",
                    mesh_.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    false
                );

                if (points0IO.headerOk())
                {
                    points0Ptr_.reset
                    (
                        new pointIOField(points0IO)
                    );
                }
                else
                {
                    points0Ptr_.reset
                    (
                        new pointIOField
                        (
                            IOobject
                            (
                                "points",
                                mesh_.facesInstance(),
                                polyMesh::meshSubDir,
                                mesh_,
                                IOobject::MUST_READ,
                                IOobject::AUTO_WRITE,
                                false
                            )
                        )
                    );
                    points0Ptr_->rename("points0");
                }
                points0Ptr_->checkIn();
                mesh_.movePoints(points0Ptr_());


                //- If the mesh is stored in constant, no motion has
                //  occurred so check if there is a version of pointD
                typeIOobject<pointVectorField> pointDIO
                (
                    pointDPtr_->name(),
                    mesh_.facesInstance(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    false
                );
                if (pointDIO.headerOk())
                {
                    pointDPtr_() ==
                        pointVectorField
                        (
                            pointDIO,
                            pointMesh::New(mesh_)
                        );
                    pointDPtr_->instance() = mesh_.time().name();
                }
            }
            pointDPtr_() == Zero;
            pointDPtr_->oldTime();
            pointDPtr_->storeOldTimes();
        }
    }

    moveMesh(FINAL_ITER);

    accelerationSchemes_.clear();

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

    if (moving_)
    {
        pointDPtr_->storeOldTimes();

        if (iter == FINAL_ITER)
        {
            points0Ptr_->write();
            pointDPtr_->write();
            mesh_.lookupObject<pointIOField>("points").write();
        }
    }
}


void Foam::regionSolvers::fluid::initialise()
{
    moveMesh(FINAL_ITER);

    // Make sure oldTime field is initialized
    if (moving_)
    {
        pointDPtr_->oldTime();
    }

    // Set old points
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
        mesh_.resetMotion();
    }
}


bool Foam::regionSolvers::fluid::changeMesh()
{
    bool changed = regionSolver::changeMesh();
    if (changed && moving_)
    {
        diffusivityPtr_->correct();

        points0Ptr_->instance() = mesh_.facesInstance();

        // // Make sure points0 are sync'd
        //TODO
        // if (Pstream::parRun())
        // {
        //     fvMeshBalance::pushUntransformedData(mesh_, points0Ptr_());
        // }

//         // Move points to the current deformed position to ensure the
//         // save old points correspond to the exiting configuration
//         tmp<pointField> tcurPoints
//         (
//             points0Ptr_() + pointDPtr_->primitiveField()
//         );
//         if (Pstream::parRun())
//         {
//             fvMeshBalance::pushUntransformedData(mesh_, tcurPoints.ref());
//         }
//         twoDPointCorrector::New(mesh_).correctPoints(tcurPoints.ref());
//         mesh_.movePoints(tcurPoints());
    }
    return changed;
}


bool Foam::regionSolvers::fluid::moveMesh(const IterType iter)
{
    regionSolver::moveMesh(iter);

    if (moving_)
    {
        pointDPtr_->oldTime();

        storePrevIter();

        // Solve point motion

        // The points have moved so before interpolation update
        // the fvMotionSolver accordingly
//         mesh_.movePoints(points0Ptr_());

        diffusivityPtr_->correct();
//         pointDPtr_->correctBoundaryConditions();
        pointDPtr_->boundaryFieldRef().updateCoeffs();

        if (iter != FINAL_ITER)
        {
            accelerationSchemes_.relax(regions_.iterNo());
        }
        else
        {
            accelerationSchemes_.updateError();
        }

        wordList patchFieldTypes
        (
            pointDPtr_->mesh().boundary().size(),
            calculatedPointPatchField<vector>::typeName
        );
        const labelList& coupledPatches = globalBoundary_.coupledPatches();
        forAll(coupledPatches, pi)
        {
            const label patchi = coupledPatches[pi];
            if
            (
                isA<fixedValuePointPatchVectorField>
                (
                    pointDPtr_->boundaryField()[patchi]
                )
            )
            {
                patchFieldTypes[patchi] =
                    valuePointPatchVectorField::typeName;
            }
        }
        //- Save boundary values
        pointVectorField::Boundary bpRelaxed
        (
            pointDPtr_->mesh().boundary(),
            pointDPtr_(),
            patchFieldTypes
        );
        bpRelaxed == pointDPtr_->boundaryField();

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

            pointVectorField::Boundary& bpointD =
                pointDPtr_->boundaryFieldRef();
            forAll(coupledPatches, pi)
            {
                const label patchi = coupledPatches[pi];
                if (isA<valuePointPatchVectorField>(bpointD[patchi]))
                {
                    valuePointPatchVectorField& ppointD =
                        dynamicCast<valuePointPatchVectorField>
                        (
                            bpointD[patchi]
                        );
                    ppointD ==
                        dynamicCast<const Field<vector>>
                        (
                            bpRelaxed[patchi]
                        );
                    ppointD.valuePointPatchVectorField::evaluate();
                }
            }
        }
//         pointDPtr_->correctBoundaryConditions();
        pointConstraints::New(pointDPtr_->mesh()).constrain
        (
            pointDPtr_(),
            false
        );

        tmp<pointField> tcurPoints
        (
            points0Ptr_() + pointDPtr_->primitiveField()
        );
        //TODO
        // if (Pstream::parRun())
        // {
        //     fvMeshBalance::pushUntransformedData(mesh_, tcurPoints.ref());
        // }
        twoDPointCorrector::New(mesh_).correctPoints(tcurPoints.ref());
        mesh_.movePoints(tcurPoints());

//         if (mesh_.moving() && (debug || regionSolver::debug))
        {
            Info<<"Mesh boundary velocity (max/mean): " << endl;
            forAll(mesh_.boundary(), patchi)
            {
                const fvPatch& p = mesh_.boundary()[patchi];
                if (!p.coupled() && returnReduce(p.size(), sumOp<label>()))
                {
                    const polyPatch& pp = p.patch();
                    const pointField& oldPoints = mesh_.oldPoints();

                    vectorField oldFc(pp.size());
                    forAll(oldFc, i)
                    {
                        oldFc[i] = pp[i].centre(oldPoints);
                    }

                    const scalar deltaT = mesh_.time().deltaTValue();

                    const vectorField Up((pp.faceCentres() - oldFc)/deltaT);

                    const volVectorField& U =
                        mesh_.lookupObject<volVectorField>("U");

                    scalarField phip
                    (
                        p.patchField<surfaceScalarField, scalar>(fvc::meshPhi(U))
                    );

                    const vectorField n(p.nf());
                    const scalarField& magSf = p.magSf();
                    scalarField Un(phip/(magSf + vSmall));
                    const vectorField pU(Up + n*(Un - (n & Up)));

                    Info<< "    " << mesh_.boundary()[patchi].name()<<": "
                        << gMaxMagSqr(pU) << "/" << gAverage(pU) << endl;
                }
            }
        }
    }

    accelerationSchemes_.print(Info);

    if (mesh_.moving())
    {
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
    return mesh_.moving();
}


// ************************************************************************* //
