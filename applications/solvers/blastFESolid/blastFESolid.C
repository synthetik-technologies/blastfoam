
#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "feMesh1.H"
#include "meshTools.H"
#include "pointFields.H"
#include "syncTools.H"
#include "DynamicList.H"
#include "materialModel.H"
#include "wedgePolyPatch.H"

#include "valuePointPatchFields.H"
#include "pointConstraints.H"
#include "twoDPointCorrector.H"

#include "displacementConstraint.H"
#include "boundaryTraction.H"

using namespace Foam;

scalar kineticEnergy(const feMesh1& mesh, const vectorField& U);

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    argList::addOption("o", "Integration order");
    argList::addOption("to", "Time integration order");
    argList::addOption("d", "Damping coefficient");
    argList::addOption("n", "Number of old KEs to save");
    argList::addBoolOption("linear", "Linear geometry");

    argList::addBoolOption
    (
        "initialize",
        "Initialize solution and write to start time"
    );
    argList::addOption("endTime", "End time of initialization");
    argList::addOption("tolerance", "Force residual");
    argList::addOption("overwrite", "Write to start time");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    feMesh1 femesh(mesh, args.optionLookupOrDefault<label>("o", 2));
    const pointField& nodes = femesh.nodes();

    scalar minDx = great;
    if (mesh.nGeometricD() == 3)
    {
        forAll(femesh.elements(), ei)
        {
            const UIndirectList<vector> pts(nodes, femesh.elements()[ei]);
            forAll(pts, i)
                for (label j = i+1; j < pts.size(); j++)
                    minDx = min(mag(pts[i] - pts[j]), minDx);
        }
    }
    else
    {
        forAll(mesh.boundaryMesh(), patchi)
        {
            if
            (
                isA<wedgePolyPatch>(mesh.boundaryMesh()[patchi])
             || isA<emptyPolyPatch>(mesh.boundaryMesh()[patchi])
            )
            {
                forAll(mesh.boundaryMesh()[patchi], ei)
                {
                    const UIndirectList<vector> pts
                    (
                        nodes,
                        femesh.boundary()[patchi].elements()[ei]
                    );
                    forAll(pts, i)
                    {
                        for (label j = i+1; j < pts.size(); j++)
                        {
                            minDx = min(mag(pts[i] - pts[j]), minDx);
                        }
                    }
                }
            }
        }
    }

    const scalar dampingCoeff = args.optionLookupOrDefault("d", -1.0);
    const label nDamp = args.optionLookupOrDefault("n", 1);

    bool initialize = args.optionFound("initialize");
    bool overwrite = initialize && args.optionFound("overwrite");
    bool write = !overwrite;
    scalar tolerance = 0.0;
    word dictName("mechanicalProperties");
    if (initialize)
    {
        dictName += ".initialize";
        if (args.optionFound("endTime"))
        {
            runTime.setEndTime(args.optionRead<scalar>("endTime"));
        }
        if (args.optionFound("tolerance"))
        {
            tolerance = args.optionRead<scalar>("tolerance");
        }
    }

    // Set time integration coefficients
    const scalar tOrder = args.optionLookupOrDefault<label>("to", 2);
    List<Pair<scalar>> tCoeffs;
    if (tOrder == 1)
    {
        tCoeffs.setSize(1);
        tCoeffs[0] = {1.0, 1.0};
    }
    else if (tOrder == 2)
    {
        tCoeffs.setSize(2);
        tCoeffs[0] = {0.0, 0.5};
        tCoeffs[1] = {1.0, 0.5};
    }
    else if (tOrder == 3)
    {
        tCoeffs.setSize(3);
        tCoeffs[0] = {1.0, -1.0/24.0};
        tCoeffs[1] = {-2.0/3.0, 3.0/4.0};
        tCoeffs[2] = {2.0/3.0, 7.0/24.0};
    }
    else if (tOrder == 4)
    {
        const scalar twoMCbrtTwo = 2.0 - Foam::cbrt(2.0);
        const scalar c14 = 0.5/twoMCbrtTwo;
        const scalar c23 = (1.0 - Foam::cbrt(2.0))/2.0/twoMCbrtTwo;
        const scalar d13 = 1.0/twoMCbrtTwo;

        tCoeffs.setSize(4);
        tCoeffs[0] = {c14, d13};
        tCoeffs[1] = {c23, -Foam::cbrt(2.0)/twoMCbrtTwo};
        tCoeffs[2] = {c23, d13};
        tCoeffs[3] = {c14, 0.0};
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported time integration order, orders 1 - 4 are "
            << "currently supported" << endl
            << abort(FatalError);
    }

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            dictName,
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    pointVectorField D
    (
        IOobject
        (
            "pointD",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        femesh.pMesh(),
        dimensionedVector(dimLength, Zero)
    );
    pointVectorField U
    (
        IOobject
        (
            "pointU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        femesh.pMesh(),
        dimensionedVector(dimVelocity, Zero)
    );
    pointVectorField force
    (
        IOobject
        (
            "force",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        femesh.pMesh(),
        dimensionedVector(dimAcceleration*dimDensity, Zero)
    );
    pointScalarField pressure
    (
        IOobject
        (
            "pressure",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        femesh.pMesh(),
        dimensionedScalar(dimPressure, Zero)
    );
    pointScalarField::Internal M
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh
        ),
        femesh.pMesh(),
        dimensionedScalar(dimDensity, Zero)
    );
    tmp<pointScalarField::Internal> M0ByM;

    pointSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        femesh.pMesh(),
        dimensionedSymmTensor(dimPressure, Zero)
    );


    PtrList<entry> materials(mechanicalProperties.lookup("mechanical"));
    const bool planeStress = mechanicalProperties.lookup<bool>("planeStress");
    materialModel::GeoType geoType =
        materialModel::TOTAL_LAGRANGIAN;
    if (args.optionFound("linear"))
    {
        geoType = materialModel::LINEAR;
    }

    autoPtr<materialModel> materialPtr
    (
        materialModel::New
        (
            materials[0].dict(),
            femesh,
            D,
            U,
            planeStress,
            geoType
        )
    );
    materialModel& material = materialPtr();

    PtrList<displacementConstraint> constraints
    (
        mechanicalProperties.lookup("displacementConstraints"),
        displacementConstraint::iNew(D, U)
    );

    PtrList<boundaryTraction> loads
    (
        mechanicalProperties.lookup("boundaryLoads"),
        boundaryTraction::iNew
        (
            femesh,
            geoType == materialModel::TOTAL_LAGRANGIAN
          ? &D
          : nullptr
        )
    );

    const scalar rho = material.rho().value();
    scalar maxVp = material.vp();

    // Mass matrix (diagonal)
    forAll(femesh.elements(), ei)
    {
        const element& elem = femesh.elements()[ei];

        UIndirectList<scalar> m_loc(M, elem);

        const integrationRule& ir = elem.ir();
        forAll(ir, rulei)
        {
            const integrationPoint& ip = ir[rulei];

            const scalarList& shape = femesh.shapes()[ei][rulei];
            const scalar w = ip.w()*femesh.Ws()[ei][rulei]*rho;
            forAll(shape, si)
            {
                forAll(shape, sj)
                {
                    m_loc[si] += w*shape[si]*shape[sj];
                }
            }
        }
    }

    syncTools::syncPointList(mesh, M, plusEqOp<scalar>(), 0.0);
    if (gMin(M) < small)
    {
        FatalErrorInFunction
            << "Mass matrix is singular" << endl
            << abort(FatalError);
    }


    scalar KE = 0.0;
    DynamicList<scalar> KEold;
    while (runTime.running())
    {
        if (runTime.controlDict().lookup<bool>("adjustTimeStep"))
        {
//             tmp<pointField> tdisplacedPoints;
//             if (geoType == materialModel::TOTAL_LAGRANGIAN)
//             {
//                 tdisplacedPoints = mesh.points() + D.primitiveField();
//             }
//             else
//             {
//                 tdisplacedPoints = tmp<pointField>(mesh.points());
//             }
//             const pointField& displacedPoints = tdisplacedPoints();
//
//             minDx = great;
//             if (mesh.nGeometricD() == mesh.nSolutionD())
//             {
//                 forAll(femesh.elements(), ei)
//                 {
//                     const UIndirectList<vector> pts
//                     (
//                         displacedPoints,
//                         femesh.elements()[ei]
//                     );
//                     forAll(pts, i)
//                     {
//                         for (label j = i+1; j < pts.size(); j++)
//                         {
//                             minDx = min(mag(pts[i] - pts[j]), minDx);
//                         }
//                     }
//                 }
//             }
//             else
//             {
//                 forAll(mesh.boundaryMesh(), patchi)
//                 {
//                     if
//                     (
//                         isA<wedgePolyPatch>(mesh.boundaryMesh()[patchi])
//                      || isA<emptyPolyPatch>(mesh.boundaryMesh()[patchi])
//                     )
//                     {
//                         forAll(mesh.boundaryMesh()[patchi], ei)
//                         {
//                             const UIndirectList<vector> pts
//                             (
//                                 displacedPoints,
//                                 femesh.boundary()[patchi].elements()[ei]
//                             );
//                             forAll(pts, i)
//                             {
//                                 for (label j = i+1; j < pts.size(); j++)
//                                 {
//                                     minDx =
//                                         min
//                                         (
//                                             mag(pts[i] - pts[j]),
//                                             minDx
//                                         );
//                                 }
//                             }
//                         }
//                     }
//                 }
//             }
//             reduce(minDx, minOp<scalar>());

            scalar CoNum = runTime.deltaTValue()*maxVp/minDx;
            const scalar maxCo =
                runTime.controlDict().lookup<scalar>("maxCo");
            double newDeltaT = maxCo*runTime.deltaTValue()/CoNum;
            runTime.setDeltaT(newDeltaT);
            Info << nl << "Max Courant Number: "<< CoNum << endl;
        }
        runTime++;
        Info<<"Time " << runTime.timeName()
            << ", deltaT = " << runTime.deltaTValue() << endl;


        maxVp = 0.0;
        forAll(tCoeffs, stepi)
        {
            Info<< "Sympletic step " << stepi << endl;

            force = Zero;
            sigma = Zero;
            pressure = Zero;

            // Add force to all elements based on material model
            material.preUpdate();
            forAll(femesh.elements(), ei)
            {
                maxVp = max(maxVp, material.addForce(force, sigma, ei));
            }
            reduce(maxVp, maxOp<scalar>());

            // Apply boundary tractions
            forAll(loads, i)
            {
                loads[i].update();
                const polyPatch& patch =
                    mesh.boundaryMesh()[loads[i].patchID()];
                const fePatch1& fepatch =
                    femesh.boundary()[patch.index()];
                forAll(patch, ei)
                {
                    //- Local nodes
                    const element& elem = fepatch.elements()[ei];

                    const UIndirectList<vector> d(D, elem);
                    UIndirectList<vector> f_loc(force, elem);
                    UIndirectList<scalar> p_loc(pressure, elem);

                    const integrationRule& ir = elem.ir();

                    forAll(ir, rulei)
                    {
                        const List<scalar>& shape =
                            fepatch.shapes()[ei][rulei];
                        const vector n
                        (
                            geoType == materialModel::TOTAL_LAGRANGIAN
                          ? normalised
                            (
                                fepatch.nodeNormals()[ei][rulei]
                              & material.calcF
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

                        const vector t
                        (
                            loads[i].traction(nodeLabels, shape, n)
                        );
                        const scalar p =
                            loads[i].pressure(nodeLabels, shape, n);
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
                mesh,
                force,
                plusEqOp<vector>(),
                vector::zero
            );
            pointConstraints::New(femesh.pMesh()).constrain(force);


            // Sum stress on coupled points
            syncTools::syncPointList
            (
                mesh,
                sigma,
                plusEqOp<symmTensor>(),
                symmTensor::zero
            );
            // Scale stress by row summed weights
            sigma /= femesh.W();

    //         // Sum stress on coupled points
    //         syncTools::syncPointList
    //         (
    //             mesh,
    //             pressure,
    //             plusEqOp<scalar>(),
    //             0.0
    //         );
    //         pressure /= femesh.W();

            material.postUpdate(tCoeffs[stepi].first());

            // Update velocity
            U.ref() +=
                tCoeffs[stepi].second()
               *force()
               *runTime.deltaT()
               /M;

            // Update displacement
            D.storeOldTimes();
            D.primitiveFieldRef() +=
                tCoeffs[stepi].first()
               *runTime.deltaTValue()
               *U.primitiveField();

            // Override velocity patches where the displacement is fixed
            forAll(constraints, i)
            {
                constraints[i].constrain();

                // Zero force on constrained nodes
                UIndirectList<vector>(force, constraints[i].nodes()) = Zero;
            }
        }

        // Damp solution
        if (dampingCoeff >= 0)
        {
            KE = kineticEnergy(femesh, U);
            Info<<"Kinetic energy: " << KE << endl;
            if (KEold.size() < nDamp)
            {
                if (!KEold.size() || KE < KEold.last())
                {
                    KEold.append(KE);
                }
                else
                {
                    KEold.clear();
                }
            }
            else
            {
                if (KE < KEold.last())
                {
                    U *= dampingCoeff;
                    KEold.clear();
                }
                else
                {
                    for (label i = 1; i < KEold.size(); i++)
                    {
                        KEold[i-1] = KEold[i];
                    }
                    KEold.last() = KE;
                }
            }
        }
        U.correctBoundaryConditions();


        if (write)
        {
            runTime.write();
        }
        else
        {
            scalar forceMaxMagSqr(mag(gMaxMagSqr(force.primitiveField())));
            Info<<"Max force: " << forceMaxMagSqr << endl;

            if (forceMaxMagSqr < tolerance)
            {
                runTime.stopAt(Time::stopAtControl::nextWrite);
            }
        }

        Info<< "Min/max |D|: " << minMagSqr(D.primitiveField())
            << ", " << maxMagSqr(D.primitiveField()) << nl << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    if (initialize)
    {
        runTime.setTime(runTime.startTime().value(), 0);
        runTime.writeNow();
    }
    Info<< nl << "done" <<endl;
    return 0;
}


scalar kineticEnergy(const feMesh1& mesh, const vectorField& U)
{
    scalar KE = 0.0;
    const elementList& elements = mesh.elements();
    const List<List<scalarList>>& shapesI = mesh.shapes();
    const List<List<scalar>>& Ws = mesh.Ws();

    forAll(elements, ei)
    {
        const element& elem = elements[ei];
        const UIndirectList<vector> u(U, elem);
        forAll(elem.ir(), rulei)
        {
            const List<scalar>& shape = shapesI[ei][rulei];
            scalar magSqrU_ip(Zero);
            forAll(shape, si)
            {
                magSqrU_ip += shape[si]*magSqr(u[si]);
            }
            KE += elem.ir()[rulei].w()*Ws[ei][rulei]*magSqrU_ip;
        }
    }
    return returnReduce(KE, sumOp<scalar>());
}
