/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | SYnthetik Applied Technologies
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

Application
    rotateFields

Description
    Rotate fields from 1-D to 2-D or 2-D to 3-D. Only for
    axisymmetric cases. Currently mapping is not conservative so only the
    nearest cell value (based on cell centres) will be used. Optionally,
    refinement can be used.

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "wedgeFvPatch.H"
#include "IOobjectList.H"
#include "HashSet.H"
#include "UautoPtr.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "genericFvPatchField.H"

#include "errorEstimator.H"
#include "polyMeshRefiner.H"
#include "blastMeshTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const wordHashSet ignoreField({"Cc", "Cf", "Sf", "magSf"});
bool parRun = false;
void setParRun(const bool par)
{
    Pstream::parRun() = par;
}
void resetParRun()
{
    Pstream::parRun() = parRun;
}

struct cellInfo
{
    label proc;
    label cellID;
    scalar weight;
};
typedef DynamicList<cellInfo> cellInfoList;

template<class Type>
wordList createBoundaryTypes
(
    const GeometricField<Type, fvPatchField, volMesh>& src,
    const fvMesh& targetMesh
)
{
    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        srcBoundary = src.boundaryField();

    HashTable<word> srcBoundaryTypes;
    forAll(srcBoundary, patchi)
    {
        if (!isA<genericFvPatchField<Type>>(srcBoundary[patchi]))
        {
            srcBoundaryTypes.insert
            (
                srcBoundary[patchi].patch().name(),
                srcBoundary[patchi].type()
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unknown patch type for patch "
                << srcBoundary[patchi].patch().name()
                << " for " << src.name() << nl
                << "Perhaps the library is missing?" << nl
                << "include the necessary library with"
                << " \'libs (\"lib*.so\")\'" << nl
                << " in the controlDict" << endl
                << abort(FatalError);
        }
    }

    wordList targetBoundaryTypes
    (
        targetMesh.boundary().size(),
        calculatedFvPatchField<Type>::typeName
    );
    forAll(targetMesh.boundary(), patchi)
    {
        if (srcBoundaryTypes.found(targetMesh.boundaryMesh()[patchi].name()))
        {
            targetBoundaryTypes[patchi] =
                srcBoundaryTypes[targetMesh.boundaryMesh()[patchi].name()];
        }
    }
    return targetBoundaryTypes;
}


template<class Type>
void mapVolFields
(
    const PtrList<fvMesh>& sourceMeshes,
    const fvMesh& targetMesh,
    const List<cellInfoList>& cellMap,
    const List<cellInfoList>& extendedCellMap,
    const IOobjectList& objects,
    const tensorField& R,
    const HashSet<word>& mapFields,
    const bool store = false
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;
    IOobjectList fields = objects.lookupClass(FieldType::typeName);
    forAllIter(IOobjectList, fields, fieldIter)
    {
        if (ignoreField.found(fieldIter()->name()))
        {
            continue;
        }
        typeIOobject<FieldType> fieldTargetIOobject
        (
            fieldIter()->name(),
            targetMesh.time().name(),
            targetMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        bool mapField = false;

        setParRun(false);
        PtrList<FieldType> fieldSources(sourceMeshes.size());
        UPtrList<const FieldType> UfieldSources(sourceMeshes.size());
        forAll(sourceMeshes, proci)
        {
            typeIOobject<FieldType> fieldSourceIOobject
            (
                fieldIter()->name(),
                sourceMeshes[proci].time().name(),
                sourceMeshes[proci],
                IOobject::MUST_READ
            );
            if (sourceMeshes[proci].foundObject<FieldType>(fieldIter()->name()))
            {
                UfieldSources.set
                (
                    proci,
                    &sourceMeshes[proci].lookupObject<FieldType>
                    (
                        fieldIter()->name()
                    )
                );
            }
            else if (fieldSourceIOobject.headerOk())
            {
                fieldSources.set
                (
                    proci,
                    new FieldType
                    (
                        fieldSourceIOobject,
                        sourceMeshes[proci]
                    )
                );
                UfieldSources.set(proci, &fieldSources[proci]);
            }
        }
        resetParRun();

        autoPtr<FieldType> fieldTargetPtr;
        UautoPtr<const List<cellInfoList>> mapPtr;
        bool exists = false;

        if (!UfieldSources.set(0))
        {
        }
        else if (targetMesh.foundObject<FieldType>(fieldIter()->name()))
        {
            mapField = true;
            fieldTargetPtr.set
            (
                &targetMesh.lookupObjectRef<FieldType>
                (
                    fieldIter()->name()
                )
            );
            mapPtr.set(&cellMap);
            exists = true;
        }
        else if (fieldTargetIOobject.headerOk())
        {
            mapField = true;
            fieldTargetPtr.set
            (
                new FieldType
                (
                    fieldTargetIOobject,
                    targetMesh
                )
            );
            mapPtr.set(&cellMap);
        }
        else if (mapFields.found(fieldTargetIOobject.name()))
        {
            mapField = true;
            fieldTargetIOobject.readOpt() = IOobject::NO_READ;
            FieldType* fieldPtr =
                new FieldType
                (
                    fieldTargetIOobject,
                    targetMesh,
                    dimensioned<Type>
                    (
                        "0",
                        UfieldSources[0].dimensions(),
                        pTraits<Type>::zero
                    ),
                    createBoundaryTypes<Type>(UfieldSources[0], targetMesh)
                );
            if (store)
            {
                fieldPtr->store(fieldPtr);
                exists = true;
            }
            fieldTargetPtr.set(fieldPtr);
            mapPtr.set(&extendedCellMap);
        }

        if (mapField)
        {
            Info<< "    mapping " << fieldIter()->name() << endl;

            // Read fieldTarget
            FieldType& fieldTarget = fieldTargetPtr();

            const List<cellInfoList>& map = mapPtr();

            forAll(map, celli)
            {
                if (map[celli].size())
                {
                    fieldTarget[celli] = Zero;
                    Type vSum = Zero;
                    scalar sumW = 0.0;
                    forAll(map[celli], j)
                    {
                        const cellInfo& info = map[celli][j];
                        vSum += info.weight*UfieldSources[info.proc][info.cellID];
                        sumW += info.weight;
                    }
                    fieldTarget[celli] = transform(R[celli], vSum/sumW);
                }
            }
            forAll(fieldTarget.boundaryField(), patchi)
            {
                fieldTarget.boundaryFieldRef()[patchi] =
                    fieldTarget.boundaryField()[patchi].patchInternalField();
            }

            if (!exists)
            {
                fieldTarget.write();
            }
            else
            {
                // Remove the pointer because it should not be deleted
                // here
                fieldTargetPtr.ptr();
            }
        }

#ifdef FULLDEBUG
        else
        {
            Info<< "    Not mapping " << fieldIter()->name() << nl
                << "         Add to \"additionalFields\" if you would "
                << "like to include it" << endl;
        }
#endif
    }
}


void mapFields
(
    const PtrList<fvMesh>& sourceMeshes,
    const fvMesh& targetMesh,
    const List<cellInfoList>& cellMap,
    const List<cellInfoList>& extendedCellMap,
    const tensorField& R,
    const HashSet<word>& additionalFields,
    const bool store = false
);


Foam::Pair<Foam::vector> calculateAxis(const fvMesh& mesh)
{
    Pair<vector> axis(vector::one, Zero);
    vector& rAxis = axis[0];
    vector& yAxis = axis[1];

    List<vector> foundAxis;
    forAll(mesh.boundaryMesh(), patchi)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchi]))
        {
            const wedgePolyPatch& wp = dynamicCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchi]
            );
            bool found = false;
            vector a = cmptMag(wp.axis());
            forAll(foundAxis, ai)
            {
                if (mag(a - foundAxis[ai]) < 1e-6)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                foundAxis.append(a);
                rAxis -= cmptMag(wp.centreNormal()) + a;
                yAxis += a;
            }
        }
    }
    forAll(yAxis, cmpti)
    {
        rAxis[cmpti] = min(rAxis[cmpti]*pos(rAxis[cmpti] - 1e-6), 1.0);
        yAxis[cmpti] = min(yAxis[cmpti]*pos(yAxis[cmpti] - 1e-6), 1.0);
    }
    return axis;
}


void calcMapAndR
(
    const PtrList<indexedOctree<treeDataCell>>& icos,
    const PtrList<fvMesh>& sourceMeshes,
    const fvMesh& targetMesh,
    const scalar& maxR,
    const vector& sourceCentre,
    const vector& targetCentre,
    const vector& rotationAxis,
    const vector& rAxis,
    const bool nearest,
    List<cellInfoList>& cellMap,
    List<cellInfoList>& extendedCellMap,
    tensorField& R
)
{
    Info<< "Calulating map and rotation tensors" << endl;
    label nSourceD = sourceMeshes[0].nGeometricD();
    vector sourceD(sourceMeshes[0].geometricD());
    vector targetD(targetMesh.geometricD());

    pointField transformedPoints(targetMesh.points());
    if (mag(rAxis) > small)
    {
        forAll(transformedPoints, pointi)
        {
            // Get the position on the target mesh
            vector ptTarget =
                cmptMultiply
                (
                    targetMesh.points()[pointi],
                    targetD
                );

            // Offset from the center of the source mesh
            vector nTarget = (ptTarget - targetCentre);

            // Radius
            scalar r = mag(nTarget);

            // Offset from the source mesh center (only solved directions)
            vector nSource(Zero);
            if (nSourceD == 1)
            {
                nSource = r*rAxis;
            }
            else
            {
                scalar y = nTarget & rotationAxis;
                scalar x = mag((nTarget - y*rotationAxis));
                nSource = y*rotationAxis + x*rAxis;
            }

            // Actual point on the source mesh
            transformedPoints[pointi] = nSource + sourceCentre;
        }
    }

    forAll(cellMap, celli)
    {
        const cell& c = targetMesh.cells()[celli];
        const vector& cc = targetMesh.cellCentres()[celli];

        // Mapping has already been set
        // if (cellMap[celli].size())
        // {
        //     continue;
        // }
        cellMap[celli].clear();
        // Get the position on the target mesh
        vector ptTarget =
            cmptMultiply
            (
                targetMesh.cellCentres()[celli],
                targetD
            );

        // Offset from the center of the source mesh
        vector nTarget = (ptTarget - targetCentre);

        // Radius
        scalar r = mag(nTarget);

        // Offset from the source mesh center (only solved directions)
        vector nSource(Zero);
        if (nSourceD == 1)
        {
            nSource = r*rAxis;
        }
        else
        {
            scalar y = nTarget & rotationAxis;
            scalar x = mag((nTarget - y*rotationAxis));
            nSource = y*rotationAxis + x*rAxis;
        }

        // Actual point on the source mesh
        vector ptSource = nSource + sourceCentre;

        // Map from the source mesh to the target mesh
        // Keep track of the maximum volume (weight) for extended cell
        scalar sumV = 0.0;
        scalar minDSqr = great;
        cellInfoList& infos = cellMap[celli];
        if (r < maxR || maxR < 0)
        {
            forAll(sourceMeshes, proci)
            {
                if (targetMesh.bounds().overlaps(sourceMeshes[proci].bounds()))
                {
                    treeBoundBox bb(c.points(targetMesh.faces(), transformedPoints));
                    forAll(sourceD, cmpti)
                    {
                        if (sourceD[cmpti] < 0)
                        {
                            bb.min()[cmpti] = -great;
                            bb.max()[cmpti] = great;
                        }
                    }

                    // Find all cells in the transformed bound box and weight based
                    // on cell volume
                    // NOTE: Ideally this would be overlap volume but that requires more
                    // work
                    labelList sCells(icos[proci].findBox(bb));
                    forAll(sCells, cj)
                    {
                        const label cellj = sCells[cj];
                        const scalar V = sourceMeshes[proci].V()[cellj];
                        const scalar dSqr =
                            magSqr(cc - sourceMeshes[proci].cellCentres()[cellj]);
                        sumV += V;
                        infos.append({proci, cellj, V});

                        if (dSqr < minDSqr)
                        {
                            minDSqr = dSqr;
                            if (infos.size() > 1)
                            {
                                Swap(infos.first(), infos.last());
                            }
                        }
                    }
                }
            }
        }

        // normalize weights
        if (nearest && infos.size())
        {
            infos.setSize(1);
            infos[0].weight = 1.0;
        }
        else if (sumV > 0)
        {
            forAll(infos, cj)
            {
                infos[cj].weight /= sumV;
            }
        }

        // Extend radius is the target point is outside of the source mesh
        scalar dist = great;
        forAll(sourceMeshes, proci)
        {
            pointIndexHit pIH = icos[proci].findNearest(ptSource, great);
            scalar curDist = mag(pIH.hitPoint() - ptSource);
            if (curDist < dist)
            {
                dist = curDist;
                extendedCellMap[celli](0) = {proci, pIH.index(), 1.0};
            }
        }
        if (r < maxR && !cellMap[celli].size())
        {
            cellMap[celli](0) = extendedCellMap[celli][0];
        }

        //- If mapping from a 2D case then there is only 1 rotation axis
        //  so to remove problems of rotating about the wrong axis we
        //  explicitly remove the axis direction from the directional
        //  vectors
        if (nSourceD == 2)
        {
            nTarget -= (nTarget & rotationAxis)*rotationAxis;
            nSource -= (nSource & rotationAxis)*rotationAxis;
        }

        // Normalise directions
        if (mag(nTarget) > small)
        {
            nTarget = nTarget/mag(nTarget);
        }
        if (mag(nSource) > small)
        {
            nSource = nSource/mag(nSource);
        }
        R[celli] = rotationTensor(nSource, nTarget);
    }
    Info<< endl;
}


//- Read and add all fields to the database
void readAllFields(const fvMesh& mesh)
{
    // Get all fields present at the current time
    // and prune the list so objects found in the registry already
    // are not double added
    IOobjectList objects;
    {
        IOobjectList tobjects(mesh, mesh.time().name());

        forAllConstIter
        (
            IOobjectList,
            tobjects,
            iter
        )
        {
            if (!mesh.foundObject<regIOobject>(iter.key()))
            {
                objects.insert
                (
                    iter.key(),
                    new IOobject(*iter())
                );
            }
        }
    }

    meshTools::readGeoFields<scalar, fvPatchField, volMesh>(mesh, objects);
    meshTools::readGeoFields<vector, fvPatchField, volMesh>(mesh, objects);
    meshTools::readGeoFields<symmTensor, fvPatchField, volMesh>(mesh, objects);
    meshTools::readGeoFields<sphericalTensor, fvPatchField, volMesh>(mesh, objects);
    meshTools::readGeoFields<tensor, fvPatchField, volMesh>(mesh, objects);

    meshTools::readGeoFields<scalar, fvsPatchField, surfaceMesh>(mesh, objects);
    meshTools::readGeoFields<vector, fvsPatchField, surfaceMesh>(mesh, objects);
    meshTools::readGeoFields<symmTensor, fvsPatchField, surfaceMesh>(mesh, objects);
    meshTools::readGeoFields<sphericalTensor, fvsPatchField, surfaceMesh>(mesh, objects);
    meshTools::readGeoFields<tensor, fvsPatchField, surfaceMesh>(mesh, objects);

    meshTools::readGeoFields<scalar, pointPatchField, pointMesh>(mesh, objects);
    meshTools::readGeoFields<vector, pointPatchField, pointMesh>(mesh, objects);
    meshTools::readGeoFields<symmTensor, pointPatchField, pointMesh>(mesh, objects);
    meshTools::readGeoFields<sphericalTensor, pointPatchField, pointMesh>(mesh, objects);
    meshTools::readGeoFields<tensor, pointPatchField, pointMesh>(mesh, objects);
}


void refine
(
    const PtrList<indexedOctree<treeDataCell>>& icos,
    const PtrList<fvMesh>& sourceMeshes,
    fvMesh& targetMesh,
    const scalar maxR,
    const vector& sourceCentre,
    const vector& targetCentre,
    const vector& rotationAxis,
    const vector& rAxis,
    const bool nearest,
    const wordList& additionalFieldNames
)
{
    dictionary refineDict;
    {
        typeIOobject<IOdictionary> dynamicMeshDictIO
        (
            IOobject
            (
                "dynamicMeshDict",
                targetMesh.time().constant(),
                targetMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        typeIOobject<IOdictionary> rotateFieldsDictIO
        (
            IOobject
            (
                "rotateFieldsDict",
                targetMesh.time().system(),
                targetMesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        if (dynamicMeshDictIO.headerOk())
        {
            IOdictionary dynamicMeshDict(dynamicMeshDictIO);
            if (dynamicMeshDict.isDict("topoChanger"))
            {
                refineDict = dynamicMeshDict.subDict("topoChanger");
            }
        }
        else
        {
            if (rotateFieldsDictIO.headerOk())
            {
                refineDict = IOdictionary(rotateFieldsDictIO);
            }
        }

        if (refineDict.empty())
        {
            WarningInFunction
                << "Refinement was specified, but neither " << nl
                << dynamicMeshDictIO.objectPath() << nl << " or " << nl
                << rotateFieldsDictIO.objectPath() << nl << " was found."
                << "Skipping." << endl;

            return;
        }
    }
    autoPtr<errorEstimator> error
    (
        errorEstimator::New(targetMesh, refineDict)
    );
    error->setForce(true);
    autoPtr<polyMeshRefiner> refiner
    (
        polyMeshRefiner::New(targetMesh, refineDict)
    );
    refiner->setForce(true);


    readAllFields(targetMesh);

    Info<< "Begining refinement iterations" << nl << endl;
    bool good = true;
    bool lastIter = false;
    label iter = 0;
    label maxIter = refineDict.lookupOrDefault("maxRotateRefineIter", 5);
    while (good)
    {
        if (++iter > maxIter)
        {
            lastIter = true;
        }

        if (lastIter)
        {
            good = false;
        }

        Info<<"Iteration " << iter << endl;

        List<cellInfoList> cellMap(targetMesh.nCells());
        List<cellInfoList> extendedCellMap(targetMesh.nCells());
        tensorField R(targetMesh.nCells(), tensor::I);

        calcMapAndR
        (
            icos,
            sourceMeshes,
            targetMesh,
            maxR,
            sourceCentre,
            targetCentre,
            rotationAxis,
            rAxis,
            nearest,
            cellMap,
            extendedCellMap,
            R
        );

        // Map fields from the source mesh to the target mesh
        mapFields
        (
            sourceMeshes,
            targetMesh,
            cellMap,
            extendedCellMap,
            R,
            additionalFieldNames,
            true
        );
        if (!lastIter)
        {
            error->update();
            lastIter =
                !refiner->refine(error->error(), error->maxRefinement());
        }

        Info<< "ExecutionTime = " << targetMesh.time().elapsedCpuTime() << " s"
            << "  ClockTime = " << targetMesh.time().elapsedClockTime() << " s"
            << nl << endl;
    }
    refiner->write();

    Info<< "Final target mesh size: " << targetMesh.nCells() << nl
        << endl;
}


void addOptions()
{
    //- Add options
    argList::addNote
    (
        "Rotational extrusion and map volume fields from one mesh to another\n"
    );
    argList::validArgs.append("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "specify the source time"
    );
    argList::addOption
    (
        "sourceRegion",
        "word",
        "specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "specify the target region"
    );
    argList::addBoolOption
    (
        "parallelSource",
        "the source is decomposed"
    );
    argList::addBoolOption
    (
        "extend",
        "Use the closest cell value if a given cell is outside of the mesh"
    );
    argList::addOption
    (
        "maxR",
        "scalar|'maximum radius'",
        "Cut off radius to map"
    );
    argList::addOption
    (
        "centre",
        "vector|'(0 0 0)'",
        "Location of center in the target mesh"
    );
    argList::addOption
    (
        "sourceCentre",
        "vector|'(0 0 0)'",
        "Location of center in the source mesh"
    );
    argList::addOption
    (
        "additionalFields",
        "wordList|'(rho U)'",
        "List of additional fields to map"
    );
    argList::addBoolOption
    (
        "uniform",
        "Copy uniform objects (not time)"
    );
    argList::addBoolOption
    (
        "refine",
        "Iteratively map, rotate and refine the mesh"
    );
    argList::addBoolOption
    (
        "tets",
        "Use cell tet decomposition"
    );
    argList::addBoolOption
    (
        "points0",
        "Write points0 field for moving meshes"
    );
    argList::addBoolOption
    (
        "nearest",
        "Use the nearest cell to map"
    );
    #include "addRegionOption.H"
}

// ************************************************************************* //
