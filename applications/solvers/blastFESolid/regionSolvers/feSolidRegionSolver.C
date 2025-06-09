/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.


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

#include "feSolidRegionSolver.H"

#include "syncTools.H"

#include "wedgePolyPatch.H"
#include "valuePointPatchFields.H"
#include "pointConstraints.H"
#include "twoDPointCorrector.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace regionSolvers
{
    defineTypeNameAndDebug(feSolid, 0);
    addToRunTimeSelectionTable(regionSolver, feSolid, dictionary);
}
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::feSolid::feSolid
(
    fvMesh& mesh,
    const regionSolverList& regions
)
:
    regionSolver(mesh, regions),

    femesh_(mesh, 2),
    solidProperties_
    (
        IOobject
        (
            "solidProperties",
            runTime_.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mechanicalProperties_
    (
        IOobject
        (
            "mechanicalProperties",
            runTime_.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    geoType_(materialModel::LINEAR),
    D_
    (
        IOobject
        (
            "pointD",
            runTime_.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        femesh_.pMesh(),
        dimensionedVector(dimLength, Zero)
    ),
    U_
    (
        IOobject
        (
            "pointU",
            runTime_.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        femesh_.pMesh(),
        dimensionedVector(dimVelocity, Zero)
    ),
    L_(femesh_.elements().size(), 0.0),
    force_
    (
        IOobject
        (
            "pointForce",
            runTime_.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        femesh_.pMesh(),
        dimensionedVector(dimAcceleration*dimDensity, Zero)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime_.name(),
            mesh
        ),
        femesh_.pMesh(),
        dimensionedScalar(dimDensity, Zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            runTime_.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        femesh_.pMesh(),
        dimensionedSymmTensor(dimPressure, Zero)
    ),
    pressure_
    (
        IOobject
        (
            "pressure",
            runTime_.name(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        femesh_.pMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    maxVp_(0.0),
    constraints_
    (
        mechanicalProperties_.lookup("displacementConstraints"),
        displacementConstraint::iNew(D_, U_)
    ),
    loads_(0)
{
    globalPolyBoundaryMesh& globalPatches =
        globalPolyBoundaryMesh::New(mesh_);
    globalPatches.setDisplacementField(mesh_.name(), "none");

    const word solidType(solidProperties_.lookup("solidModel"));
    if (solidType == "linear")
    {
        geoType_ = materialModel::LINEAR;
        globalPolyBoundaryMesh::New(mesh).setDisplacementField
        (
            mesh_.name(),
            "none"
        );
    }
    else if (solidType == "updatedLagrangian")
    {
        geoType_ = materialModel::UPDATE_LAGRANGIAN;
        globalPolyBoundaryMesh::New(mesh).setDisplacementField
        (
            mesh_.name(),
            "none"
        );
    }
    else if (solidType == "totalLagrangian")
    {
        geoType_ = materialModel::TOTAL_LAGRANGIAN;
        globalPolyBoundaryMesh::New(mesh).setDisplacementField
        (
            mesh_.name(),
            D_.name()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Invalid feSolid type " << solidType << ". "
            << "Valid feSolids are" << nl
            << "    linear" << nl
            << "    updatedLagrangian" << nl
            << "    totalLagrangian" << nl
            << abort(FatalError);
    }

    PtrList<entry> materialDicts
    (
        mechanicalProperties_.lookup("mechanical")
    );
    const bool planeStress =
        mechanicalProperties_.lookup<bool>("planeStress");
    materials_.setSize(materialDicts.size());
    forAll(materialDicts, i)
    {
        materials_.set
        (
            i,
            (
                materialModel::New
                (
                    materialDicts[i].dict(),
                    femesh_,
                    D_,
                    U_,
                    planeStress,
                    geoType_
                )
            )
        );
    }

    // Add Displacement field to track error
//     accelerationSchemes_.addField(D_);

    // Mass matrix (diagonal)
    if (materials_.size() == 1)
    {
        const materialModel& material = materials_[0];
        forAll(femesh_.elements(), elemi)
        {
            const element& elem = femesh_.elements()[elemi];

            UIndirectList<scalar> m_loc(M_, elem);

            const integrationRule& ir = elem.ir();
            forAll(ir, rulei)
            {
                const integrationPoint& ip = ir[rulei];

                const scalarList& shape = femesh_.shapes()[elemi][rulei];
                const scalar w =
                    ip.w()
                   *femesh_.Ws()[elemi][rulei]
                   *material.rho().value();
                forAll(shape, si)
                {
                    forAll(shape, sj)
                    {
                        m_loc[si] += w*shape[si]*shape[sj];
                    }
                }
            }
        }
    }
    else
    {
        forAll(materials_, mati)
        {
            materialModel& material = materials_[mati];
            const labelList& matElements = material.elements();
            forAll(matElements, ei)
            {
                const label elemi = matElements[ei];
                const element& elem = femesh_.elements()[elemi];

                UIndirectList<scalar> m_loc(M_, elem);

                const integrationRule& ir = elem.ir();
                forAll(ir, rulei)
                {
                    const integrationPoint& ip = ir[rulei];

                    const scalarList& shape =
                        femesh_.shapes()[elemi][rulei];
                    const scalar w =
                        ip.w()
                       *femesh_.Ws()[elemi][rulei]
                       *material.rho().value();
                    forAll(shape, si)
                    {
                        forAll(shape, sj)
                        {
                            m_loc[si] += w*shape[si]*shape[sj];
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList(mesh, M_, plusEqOp<scalar>(), 0.0);
    if (gMin(M_) < small)
    {
        FatalErrorInFunction
            << "Mass matrix is singular" << endl
            << abort(FatalError);
    }

    loads_ = PtrList<boundaryTraction>
    (
        mechanicalProperties_.lookup("boundaryLoads"),
        boundaryTraction::iNew
        (
            femesh_,
            geoType_ == materialModel::TOTAL_LAGRANGIAN
          ? &D_
          : nullptr
        )
    );

    forAll(femesh_.elements(), ei)
    {
        const UIndirectList<vector> pts(femesh_.nodes(), femesh_.elements()[ei]);
        scalar  minDxi = great;
        forAll(pts, i)
        {
            for (label j = i+1; j < pts.size(); j++)
            {
                minDxi = min(mag(pts[i] - pts[j]), minDxi);
            }
        }
        L_[ei] = minDxi;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionSolvers::feSolid::~feSolid()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::feSolid::initialiseMesh(const IterType)
{}


void Foam::regionSolvers::feSolid::initialiseFields()
{}


void Foam::regionSolvers::feSolid::initialise()
{
    globalPolyBoundaryMesh& globalPatches =
        globalPolyBoundaryMesh::New(mesh_);

    if (geoType_ == materialModel::TOTAL_LAGRANGIAN)
    {
        globalPatches.setDisplacementField(mesh_.name(), "pointD");
    }
    else
    {
        globalPatches.setDisplacementField(mesh_.name(), "none");
    }
    globalPatches.setInverseDisplacement(this->mesh().name(), false);
    globalPatches.clearOut();
}


bool Foam::regionSolvers::feSolid::moveMesh(const IterType iter)
{
    regionSolver::moveMesh(iter);
    return max(mag(U_.primitiveField())) > small;
}


void Foam::regionSolvers::feSolid::solve()
{
    force_ = Zero;
    sigma_ = Zero;
    maxVp_ = 0.0;

    // Add force to all elements based on material model
    if (materials_.size() == 1)
    {
        materialModel& material = materials_[0];
        material.preUpdate();
        forAll(femesh_.elements(), ei)
        {
            maxVp_ = max(maxVp_, material.addForce(force_, sigma_, ei, L_[ei]));
        }
        material.postUpdate(1.0);
    }
    else
    {
        forAll(materials_, mati)
        {
            materialModel& material = materials_[mati];
            material.preUpdate();
            const labelList& matElements = material.elements();
            forAll(matElements, ei)
            {
                const label elemi = matElements[ei];
                maxVp_ = max
                (
                    maxVp_,
                    material.addForce(force_, sigma_, elemi, L_[ei])
                );
            }
            material.postUpdate(1.0);
        }
    }

    reduce(maxVp_, maxOp<scalar>());

    // Apply boundary tractions
    forAll(loads_, i)
    {
        loads_[i].update();
        const polyPatch& patch =
            mesh_.boundaryMesh()[loads_[i].patchID()];
        const fePatch1& fepatch =
            femesh_.boundary()[patch.index()];
        forAll(patch, ei)
        {
            //- Local nodes
            const element& elem = fepatch.elements()[ei];

            const UIndirectList<vector> d(D_, elem);
            UIndirectList<vector> f_loc(force_, elem);
            UIndirectList<scalar> p_loc(pressure_, elem);

            const integrationRule& ir = elem.ir();

            forAll(ir, rulei)
            {
                const List<scalar>& shape =
                    fepatch.shapes()[ei][rulei];
                const vector n
                (
                    geoType_ == materialModel::TOTAL_LAGRANGIAN
                  ? normalised
                    (
                        fepatch.nodeNormals()[ei][rulei]
                       & materialModel::calcF
                        (
                            d,
                            fepatch.dshapes()[ei][rulei],
                            fepatch.invJs()[ei][rulei]
                        ).inv().T()
                    )
                  : fepatch.nodeNormals()[ei][rulei]
                );

                const labelList& nodeLabels =
                    fepatch.localElements()[ei];

                const scalar w = ir[rulei].w()*fepatch.Ws()[ei][rulei];

                const vector t(loads_[i].traction(nodeLabels, shape, n));
                const scalar p = loads_[i].pressure(nodeLabels, shape, n);
                const vector f(t - n*p);
                forAll(shape, si)
                {
                    f_loc[si] += f*shape[si]*w;
                    p_loc[si] = p;//p*shape[si];
                }
            }
        }
    }

    // Sum force on coupled points
    syncTools::syncPointList
    (
        mesh_,
        force_,
        plusEqOp<vector>(),
        vector::zero
    );

    // Constrain force on wegde points
    forAll(mesh_.boundaryMesh(), patchi)
    {
        if (isA<wedgePolyPatch>(mesh_.boundaryMesh()[patchi]))
        {
            const wedgePolyPatch& patch =
                dynamicCast<const wedgePolyPatch>
                (
                    mesh_.boundaryMesh()[patchi]
                );
            const vector& n = patch.centreNormal();
            const tensor& R = patch.faceT();
            const tensor invR = R.inv();
            const labelList& meshPoints = patch.meshPoints();
            forAll(meshPoints, pi)
            {
                const label pointi = meshPoints[pi];

                // Inverse rotation from wegde plane to centre plane
                vector f = transform(invR, force_[pointi]);

                // Force normal to the centre plane
                const vector fn((f & n)*n);

                // Transform force back to the wegde plane
                force_[pointi] = transform(R, f - fn);
            }
        }
    }

    // Sum stress on coupled points
    syncTools::syncPointList
    (
        mesh_,
        sigma_,
        plusEqOp<symmTensor>(),
        symmTensor::zero
    );

    // Scale stress by row summed weights
    sigma_ /= femesh_.W();


    // Update velocity
    if (runTime_.timeIndex() == 1)
    {
        U_.internalFieldRef() =
            U_.oldTime()() + force_()*runTime_.deltaT()*0.5/M_;
    }
    else
    {
        U_.internalFieldRef() =
            U_.oldTime()() + force_()*runTime_.deltaT()/M_;
    }

    // Make sure no movement in empty directions
    forAll(mesh_.solutionD(), cmpti)
    {
        if (mesh_.solutionD()[cmpti] < 0)
        {
            U_.primitiveFieldRef().replace(cmpti, 0.0);
        }
    }
//     const pointConstraints& pc =  pointConstraints::New(femesh_.pMesh());

//     pc.constrain(U_);
//     Field<vector> DD(U_.primitiveField()*runTime_.deltaTValue());
//     twoDPointCorrector::New(mesh_).correctDisplacement
//     (
//         (mesh_.points() + D_.primitiveField())(),
//         DD
//     );
//     U_.primitiveFieldRef() = DD/runTime_.deltaTValue();

    // Update displacement
    D_.primitiveFieldRef() =
        D_.oldTime().primitiveField()
      + U_.primitiveField()*runTime_.deltaTValue();
//     pc.constrainDisplacement(D_);

     // Override velocity patches where the displacement is fixed
    forAll(constraints_, i)
    {
        constraints_[i].constrain();
    }

    accelerationSchemes_.updateError();

//     const volVectorField& D = solid_->solutionD();
//
//     vector forceSum = Zero;
//     Info<< "External forces:" << incrIndent << endl;
//     forAll(D.boundaryField(), patchi)
//     {
//         const fvPatchVectorField& pD = D.boundaryField()[patchi];
//         if (isA<coupledSolidTractionFvPatchVectorField>(pD))
//         {
//             const coupledSolidTractionFvPatchVectorField& cst =
//                 dynamicCast<const coupledSolidTractionFvPatchVectorField>(pD);
//             forceSum += cst.force();
//             Info<< indent << pD.patch().name() << ":" << nl << incrIndent
//                 << indent << "solid = " << cst.force() << nl
//                 << indent << "fluid = " << cst.forceNbr() << decrIndent << endl;
//         }
//         else if (isA<solidTractionFvPatchVectorField>(pD))
//         {
//             const solidTractionFvPatchVectorField& st =
//                 dynamicCast<const solidTractionFvPatchVectorField>(pD);
//             forceSum += st.force();
//             Info<< indent << pD.patch().name() << ": "
//                 << st.force() << endl;
//         }
//     }
//     Info<< indent << "Total: " << forceSum << decrIndent << nl << endl;
//
//     // Turn solver information back on
//     SolverPerformance<vector>::debug = 1;
}


void Foam::regionSolvers::feSolid::clear(const bool full)
{
    regionSolver::clear(full);
    if (geoType_ == materialModel::UPDATE_LAGRANGIAN)
    {
        mesh_.movePoints
        (
            mesh_.points()
          + U_.primitiveField()*runTime_.deltaTValue()
        );

        // Mass matrix (diagonal)
        M_ = Zero;
        if (materials_.size() == 1)
        {
            const materialModel& material = materials_[0];
            forAll(femesh_.elements(), elemi)
            {
                const element& elem = femesh_.elements()[elemi];

                UIndirectList<scalar> m_loc(M_, elem);

                const integrationRule& ir = elem.ir();
                forAll(ir, rulei)
                {
                    const integrationPoint& ip = ir[rulei];

                    const scalarList& shape =
                        femesh_.shapes()[elemi][rulei];
                    const scalar w =
                        ip.w()
                       *femesh_.Ws()[elemi][rulei]
                       *material.rho().value();
                    forAll(shape, si)
                    {
                        forAll(shape, sj)
                        {
                            m_loc[si] += w*shape[si]*shape[sj];
                        }
                    }
                }
            }
        }
        else
        {
            forAll(materials_, mati)
            {
                materialModel& material = materials_[mati];
                const labelList& matElements = material.elements();
                forAll(matElements, ei)
                {
                    const label elemi = matElements[ei];
                    const element& elem = femesh_.elements()[elemi];

                    UIndirectList<scalar> m_loc(M_, elem);

                    const integrationRule& ir = elem.ir();
                    forAll(ir, rulei)
                    {
                        const integrationPoint& ip = ir[rulei];

                        const scalarList& shape =
                            femesh_.shapes()[elemi][rulei];
                        const scalar w =
                            ip.w()
                           *femesh_.Ws()[elemi][rulei]
                           *material.rho().value();
                        forAll(shape, si)
                        {
                            forAll(shape, sj)
                            {
                                m_loc[si] += w*shape[si]*shape[sj];
                            }
                        }
                    }
                }
            }
        }

        syncTools::syncPointList(mesh_, M_, plusEqOp<scalar>(), 0.0);
        if (gMin(M_) < small)
        {
            FatalErrorInFunction
                << "Mass matrix is singular" << endl
                << abort(FatalError);
        }
    }
}


Foam::scalar Foam::regionSolvers::feSolid::CoNum() const
{
    const polyMesh& pMesh = mesh();
    tmp<pointField> tdisplacedPoints;
    if (geoType_ == materialModel::TOTAL_LAGRANGIAN)
    {
        tdisplacedPoints = mesh_.points() + D_.primitiveField();
    }
    else
    {
        tdisplacedPoints = tmp<pointField>(mesh_.points());
    }
    const pointField& displacedPoints = tdisplacedPoints();

    scalar minDx = great;

    if (pMesh.nGeometricD() == pMesh.nSolutionD())
    {
        forAll(pMesh.cells(), celli)
        {
            const cell& c = pMesh.cells()[celli];
            const scalar V = c.mag
            (
                displacedPoints,
                pMesh.faces()
            );
            scalar maxA = 0.0;
            forAll(c, fi)
            {
                maxA = max(maxA, pMesh.faces()[c[fi]].mag(displacedPoints));
            }
            minDx = min(minDx, V/maxA);
        }
    }
    else
    {
        label patchID = -1;
        forAll(pMesh.boundaryMesh(), patchi)
        {
            if (isA<wedgePolyPatch>(pMesh.boundaryMesh()[patchi]))
            {
                patchID = patchi;
                break;
            }
        }

        const polyPatch& patch = pMesh.boundaryMesh()[patchID];
        forAll(patch, fi)
        {
            const label facei = patch.start() + fi;
            const face& f = pMesh.faces()[facei];
            const labelList& fe = pMesh.faceEdges()[facei];
            const scalar A = f.mag(displacedPoints);
            scalar maxL = 0.0;
            forAll(fe, ei)
            {
                maxL =
                    max(maxL, pMesh.edges()[fe[ei]].mag(displacedPoints));
            }
            minDx = min(minDx, A/maxL);
        }
    }
    reduce(minDx, minOp<scalar>());
    scalar co = runTime_.deltaTValue()*maxVp_/minDx;

    Info<< "Courant Number ";
    if (mesh().name() != polyMesh::defaultRegion)
    {
        Info<< "for region " << mesh().name() << " ";
    }
    Info<< "Max = " << co << endl;

    return co;
}


Foam::scalar Foam::regionSolvers::feSolid::maxCo() const
{
    return
        runTime_.controlDict().lookupOrDefault
        (
            mesh_.name() + "MaxCo",
            runTime_.controlDict().lookup<scalar>("maxCo")
        );
}


Foam::scalar Foam::regionSolvers::feSolid::newDeltaT() const
{
    return regionSolver::newDeltaT();
}


// ************************************************************************* //
