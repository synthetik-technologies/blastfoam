/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 Synthetik Applied Technologies
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

Description:
    Utility to create VTKs from given patches. Optional interpolation between
    time steps

Usage:
    blastToVTK ground -fields '(p U)' -interpolationScheme cubicClamp -dt 1.1


\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeSelector.H"
#include "lookupTables1D.H"
#include "vtkWriteOps.H"
#include "PatchToPatchInterpolation.H"
#include "standAlonePatch.H"
#include "vtkTimeSeries.H"
#include "IFstream.H"
#include "PatchTools.H"
#include "IOobjectList.H"

#include "globalPolyPatch.H"
#include "amiPatchToPatchMapping.H"
#include "primitivePatchToPatchMapping.H"

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"

using namespace Foam;
// typedef PatchToPatchInterpolation<standAlonePatch, standAlonePatch>
    // Interpolation;
// typedef AMIInterpolation Interpolation;
typedef patchToPatchMapping Interpolation;

template<class Type>
void writeField
(
    std::ostream& os,
    const bool binary,
    const word& fieldName,
    const UList<Type>& fld
)
{
     const label nCmpts = pTraits<Type>::nComponents;

    os  << fieldName << ' ' << nCmpts
        << ' ' << fld.size() << ' '
        << "float" << nl;
    label fi = 0;
    List<float> data(nCmpts*fld.size());
    forAll(fld, i)
    {
        for (direction cmpti = 0; cmpti < nCmpts; cmpti++)
        {
            data[fi++] = component(fld[i], cmpti);
        }
    }

    vtkWriteOps::write(os, binary, data);
}

template<class Mesh>
const Mesh& getMesh(const fvMesh& mesh)
{
    return mesh;
}
template<>
const pointMesh& getMesh(const fvMesh& mesh)
{
    return pointMesh::New(mesh);
}

template<class GeoField>
bool foundGeoField
(
    const fvMesh& mesh,
    const word& fieldName
)
{
    bool good = true;
    IOobject io
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );
    fileHandler().readHeader
    (
        io,
        io.objectPath(),
        GeoField::typeName
    );
    good = io.headerClassName() == GeoField::typeName;

    if (!returnReduce(good, andOp<bool>()))
    {
        return false;
    }
    return true;
}

template<class Type, template<class> class GeoField>
bool writeGeoField
(
    std::ostream& os,
    const bool binary,
    const PtrList<fvMesh>& meshes,
    const label size,
    const word& fieldName,
    const label patchID,
    const scalarList& weights,
    const PtrList<Interpolation>& interps,
    const PtrList<globalPolyPatch>& globalPatches
)
{
    if (!interps.size())
    {
        return false;
    }

    Field<Type> pFld(size, Zero);
    forAll(interps, i)
    {
        // Check if the field has been registered
        const fvMesh& mesh = meshes[i];

        if (!foundGeoField<GeoField<Type>>(mesh, fieldName))
        {
            return false;
        }

        GeoField<Type> fld
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            getMesh<typename GeoField<Type>::Mesh>(meshes[i])
        );

        // Reduce the field to the global patch
        Field<Type> piField
        (
            globalPatches[i].patchFaceToGlobal(fld.boundaryField()[patchID])
        );

        // If not the master surface, interpolate to the master patch
        if (interps.set(i))
        {
            // pFld += interps[i].faceInterpolate(weights[i]*piField);
            pFld += interps[i].transferFaces
            (
                globalPatches[i].globalPatch(),
                weights[i]*piField
            );
        }
        else
        {
            pFld += weights[i]*piField;
        }

    }

    // Write
    if (Pstream::master())
    {
        pFld /= sum(weights);
        writeField
        (
            os,
            binary,
            fieldName,
            pFld
        );
    }
    return true;
}

template<class Type>
bool writePointField
(
    std::ostream& os,
    const bool binary,
    const PtrList<fvMesh>& meshes,
    const label size,
    const word& fieldName,
    const label patchID,
    const scalarList& weights,
    const PtrList<Interpolation>& interps,
    const PtrList<globalPolyPatch>& globalPatches,
    const List<labelList>& pointMaps
)
{
    if (!interps.size())
    {
        return false;
    }

    Field<Type> pFld(size, Zero);
    forAll(interps, i)
    {
        // Check if the field has been registered
        const fvMesh& mesh = meshes[i];
        if (!foundGeoField<PointField<Type>>(mesh, fieldName))
        {
            return false;
        }

        PointField<Type> fld
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh,
                IOobject::MUST_READ
            ),
            getMesh<typename PointField<Type>::Mesh>(meshes[i])
        );

        // Reduce the field to the master processor
        Field<Type> piField
        (
            globalPatches[i].patchPointToGlobal
            (
                fld.boundaryField()[patchID].patchInternalField()
            )
        );

        // Field<Type> piField(fld.boundaryField()[patchID].patchInternalField());
        // if (Pstream::parRun())
        // {
        //     // Collect values from all processors
        //     List<Field<Type>> gatheredValues(Pstream::nProcs());
        //     gatheredValues[Pstream::myProcNo()] = piField;
        //     Pstream::gatherList(gatheredValues);

        //     if (Pstream::master())
        //     {
        //         // Combine values into single field
        //         Field<Type> allValues
        //         (
        //             ListListOps::combine<Field<Type>>
        //             (
        //                 gatheredValues,
        //                 accessOp<Field<Type>>()
        //             )
        //         );

        //         // Renumber (point data) to correspond to merged points
        //         inplaceReorder(pointMaps[i], allValues);
        //         allValues.setSize(patches[i].globalPatch().nPoints());
        //         piField.transfer(allValues);
        //     }
        //     else
        //     {
        //         continue;
        //     }
        // }
        if (interps.set(i))
        {
            // pFld += interps[i].pointInterpolate(weights[i]*piField);
            pFld += interps[i].transferPoints
            (
                globalPatches[i].globalPatch(),
                weights[i]*piField
            );
        }
        else
        {
            pFld += weights[i]*piField;
        }
    }

    if (Pstream::master())
    {
        pFld /= sum(weights);
        writeField
        (
            os,
            binary,
            fieldName,
            pFld
        );
    }
    return true;
}

template<class PointField, class FaceList>
void writePatch
(
    std::ostream& os,
    const bool binary,
    const PointField& points,
    const FaceList& faces
)
{
    // Write the points
    {
        os  << "POINTS " << points.size() << " float" << nl;
        List<floatScalar> coordinates(points.size()*3);
        forAll(points, pointi)
        {
            const point& p = points[pointi];
            forAll(p, i)
            {
                coordinates[3*pointi + i] = float(p[i]);
            }
        }
        vtkWriteOps::write(os, binary, coordinates);
    }

    // Write the faces
    if (faces.size())
    {
        label nFaceNodes = 0;
        forAll(faces, facei)
        {
            nFaceNodes += faces[facei].size();
        }
        os  << "POLYGONS " << faces.size() << ' '
            << faces.size() + nFaceNodes << nl;
        labelList data(faces.size() + nFaceNodes);
        label i = 0;
        forAll(faces, facei)
        {
            data[i ++] = faces[facei].size();
            forAll(faces[facei], facePointi)
            {
                data[i ++] = faces[facei][facePointi];
            }
        }
        vtkWriteOps::write(os, binary, data);
    }
}

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Creates a vtk.series file from the saved field data.\n"
        "Optionally a different time step size can be used and\n "
        "where intermediate steps are calculated using interpolation\n\n"
    );

    argList::validArgs.append("patch");
    timeSelector::addOptions(true, true);

    argList::addBoolOption
    (
        "ascii",
        "Write VTK in ASCII"
    );
    argList::addOption
    (
        "dt",
        "Time spacing"
    );
    argList::addOption
    (
        "interpolationScheme",
        "Scheme for interpolation"
    );
    argList::addOption
    (
        "fields",
        "Fields to sample"
    );
    argList::addBoolOption
    (
        "primitive",
        "Use primitive patch interpolation"
    );
    argList::addOption
    (
        "AMI",
        "AMI interpolation method"
    );
    argList::addOption
    (
        "outputDir",
        "Output directory"
    );

    #include "addDictOption.H"
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const bool binary = !args.optionFound("ascii");

    const word regionName =
        args.optionLookupOrDefault("region", polyMesh::defaultRegion);
    word patchName(args.argRead<word>(1));
    Info<< "Creating samples on " << string(patchName) << endl;

    instantList timeDirs = timeSelector::select0
    (
        runTime,
        args
    );

    scalarList times(timeDirs.size());
    forAll(timeDirs, ti)
    {
        times[ti] = timeDirs[ti].value();
    }

    if (!times.size())
    {
        WarningInFunction
            << "No time were selected. Exiting" << endl;
        return 0;
    }

    lookupTable1D<scalar> table
    (
        times,
        "none",
        times.size() > 1
      ? args.optionLookupOrDefault<word>("interpolationScheme", "linearClamp")
      : "floor"
    );
    table.update(runTime.value());

    scalarList sampleTimes(times);
    if (args.optionFound("dt"))
    {
        scalar tMin = sampleTimes[0];
        scalar tMax = sampleTimes.last();
        scalar dt = args.optionRead<scalar>("dt");

        label nt = (tMax - tMin)/dt + 1;
        sampleTimes.setSize(nt);
        forAll(sampleTimes, ti)
        {
            sampleTimes[ti] = tMin + dt*scalar(ti);
        }

        Info<< "Sampling from " << tMin << " to " << tMax
            << ", deltaT: " << dt << nl << endl;
    }

    PtrList<Time> runTimes(1);
    runTimes.set(0, new Time(Time::controlDictName, args));

    PtrList<fvMesh> meshes(1);
    meshes.set
    (
        0,
        new fvMesh
        (
            IOobject
            (
                regionName,
                runTimes[0].timeName(),
                runTimes[0],
                IOobject::MUST_READ
            )
        )
    );

    wordList fieldNames;
    if (args.optionFound("fields"))
    {
        fieldNames = args.optionRead<wordList>("fields");
    }
    else
    {
        // Get list of objects from processor0 database
        IOobjectList objects
        (
            meshes[0],
            timeDirs[0].name()
        );
        forAllConstIter(IOobjectList, objects, iter)
        {
            if
            (
                iter()->headerClassName() == volScalarField::typeName
             || iter()->headerClassName() == volVectorField::typeName
             || iter()->headerClassName() == volSymmTensorField::typeName
             || iter()->headerClassName() == volSphericalTensorField::typeName
             || iter()->headerClassName() == volTensorField::typeName

             || iter()->headerClassName() == surfaceScalarField::typeName
             || iter()->headerClassName() == surfaceVectorField::typeName
             || iter()->headerClassName() == surfaceSymmTensorField::typeName
             || iter()->headerClassName() == surfaceSphericalTensorField::typeName
             || iter()->headerClassName() == surfaceTensorField::typeName

             || iter()->headerClassName() == pointScalarField::typeName
             || iter()->headerClassName() == pointVectorField::typeName
             || iter()->headerClassName() == pointSymmTensorField::typeName
             || iter()->headerClassName() == pointSphericalTensorField::typeName
             || iter()->headerClassName() == pointTensorField::typeName
            )
            {
                fieldNames.append(iter.key());
            }
        }
    }
    Info<< "Sampling fields:" << endl << incrIndent;
    forAll(fieldNames, fieldi)
    {
        Info<< indent << fieldNames[fieldi] << endl;
    }
    Info<< decrIndent << endl;


    label patchID = meshes[0].boundaryMesh().findPatchID(patchName);
    if (patchID < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName
            << " not found." << nl
            << "Available patch names: "
            << meshes[0].boundaryMesh().names() << endl
            << abort(FatalError);
    }

    // Set the AMI dict
    dictionary interpDict;
    interpDict.add
    (
        "interpolationMethod",
        args.optionLookupOrDefault<word>("AMI", "faceAreaWeight")
    );

    fileName outputDir(args.optionLookupOrDefault<fileName>("outputDir", "VTK"));
    vtkTimeSeries timeSeries(outputDir, 0, true);
    forAll(sampleTimes, ti)
    {

        const scalar t = sampleTimes[ti];
        const word tName = Time::timeName(t);
        timeSeries.insert(t);
        table.update(t);

        Info<< "Sampling time = " << tName << endl;

        const labelList& indices = table.indices();
        const scalarList& weights = table.weights();

        DynamicList<label> Is(indices.size());
        DynamicList<scalar> ws(weights.size());

        if (runTimes.size() < indices.size())
        {
            runTimes.setSize(indices.size());
            meshes.setSize(indices.size());
        }

        // Track the master ID based on the highest weight
        label masterID = -1;
        scalar minDt = great;

        // Counter for actual number of indicies added
        label I = 0;
        forAll(weights, i)
        {
            if (mag(weights[i]) > 1e-10)
            {
                Is.append(indices[i]);
                ws.append(weights[i]);

                scalar dt = mag(t - times[Is[i]]);

                // Use the exact time
                if (dt < vSmall)
                {
                    masterID = I;
                    minDt = -great;
                }
                else if (dt <= minDt)
                {
                    masterID = I;
                    minDt = dt;
                }
                I++;
            }
        }

        incrIndent(Info);
        forAll(Is, i)
        {
            Info<< indent << "Time: " << timeDirs[Is[i]].value()
                << ", weight: "<< ws[i] << endl;
        }
        Info<< decrIndent << endl;

        UPtrList<const polyPatch> meshPatches(Is.size());
        forAll(meshPatches, tj)
        {
            if (!runTimes.set(tj))
            {
                runTimes.set
                (
                    tj,
                    new Time
                    (
                        Time::controlDictName,
                        args
                    )
                );
            }
            runTimes[tj].setTime(timeDirs[Is[tj]].value(), ti);

            if (!meshes.set(tj))
            {
                meshes.set
                (
                    tj,
                    new fvMesh
                    (
                        IOobject
                        (
                            regionName,
                            runTimes[tj].timeName(),
                            runTimes[tj],
                            IOobject::MUST_READ
                        )
                    )
                );
                if (patchID < 0)
                {
                    patchID = meshes[tj].boundaryMesh().findPatchID(patchName);
                    if (patchID < 0)
                    {
                        FatalErrorInFunction
                            << "Patch named " << patchName
                            << " not found." << nl
                            << "Available patch names: "
                            << meshes[0].boundaryMesh().names() << endl
                            << abort(FatalError);
                    }
                }
            }
            else
            {
                meshes[tj].readUpdate();
            }
            meshPatches.set(tj, &meshes[tj].boundaryMesh()[patchID]);
        }

        PtrList<globalPolyPatch> globalPatches(Is.size());
        List<labelList> pointMaps(Is.size());
        forAll(globalPatches, patchi)
        {
            globalPatches.set
            (
                patchi,
                new globalPolyPatch(meshPatches[patchi])
            );
        }

        const globalPolyPatch& masterPatch = globalPatches[masterID];
        const standAlonePatch& masterSurface = masterPatch.globalPatch();
        label nFaces = masterSurface.size();
        label nPoints = masterSurface.nPoints();
        PtrList<Interpolation> interps(Is.size());
        forAll(interps, patchi)
        {
            if (patchi != masterID)
            {
                // Flip the non master patches so patches are correctly orriented
                standAlonePatch& gpp =
                    const_cast<standAlonePatch&>(globalPatches[patchi].globalPatch());
                forAll(gpp, facei)
                {
                    gpp[facei].flip();
                };

                if (args.optionFound("primitive"))
                {
                    interps.set
                    (
                        patchi,

                        new patchToPatchMappings::primitivePatchToPatchMapping
                        (
                            interpDict,
                            meshPatches[patchi],
                            meshPatches[masterID],
                            globalPatches[patchi],
                            globalPatches[masterID]
                        )
                    );
                }
                else
                {
                    interps.set
                    (
                        patchi,

                        new patchToPatchMappings::amiPatchToPatchMapping
                        (
                            interpDict,
                            meshPatches[patchi],
                            meshPatches[masterID],
                            globalPatches[patchi],
                            globalPatches[masterID]
                        )
                    );
                }
            }
        }

        // Make the directory
        fileName dir = outputDir / tName;
        mkDir(dir);
        fileName file = dir / patchName + ".vtk";

        // Open the file
        std::ofstream os(file, std::ios::binary);

        // Write the header
        if (Pstream::master())
        {
            vtkWriteOps::writeHeader(os, binary, patchName);
            os << "DATASET POLYDATA" << nl;

            writePatch
            (
                os,
                binary,
                masterSurface.points(),
                masterSurface
            );
        }

        #define FoundGeoField(Type, GeoField)   \
        forAll(globalPatches, i)                \
        {                                       \
            found =                             \
                found                           \
             || foundGeoField<GeoField<Type>>   \
                (                               \
                    meshes[i],                  \
                    fieldNames[fieldi]          \
                );                              \
        }

        #define WriteGeoField(Type, GeoField)   \
        writeGeoField<Type, GeoField>           \
        (                                       \
            os,                                 \
            binary,                             \
            meshes,                             \
            nFaces,                             \
            fieldName,                          \
            patchID,                            \
            ws,                                 \
            interps,                            \
            globalPatches                       \
        );
        #define WritePointField(Type, null)     \
        writePointField<Type>                   \
        (                                       \
            os,                                 \
            binary,                             \
            meshes,                             \
            nPoints,                            \
            fieldName,                          \
            patchID,                            \
            ws,                                 \
            interps,                            \
            globalPatches,                      \
            pointMaps                           \
        );


        // Loop through all the fields and check if exists and
        // if it is a vol/surface field or a point field
        DynamicList<word> faceFields;
        DynamicList<word> pointFields;
        forAll(fieldNames, fieldi)
        {
            bool found = false;
            FOR_ALL_FIELD_TYPES(FoundGeoField, VolField);
            FOR_ALL_FIELD_TYPES(FoundGeoField, SurfaceField);
            if (found)
            {
                faceFields.append(fieldNames[fieldi]);
                found = true;
            }
            bool foundFace = found;
            found = false;
            FOR_ALL_FIELD_TYPES(FoundGeoField, PointField);
            if (found)
            {
                pointFields.append(fieldNames[fieldi]);
                found = true;
            }

            if (!foundFace && !found)
            {
                WarningInFunction
                    << "Did not find " << fieldNames[fieldi]
                    << " for time " << runTimes[masterID].timeName() << endl;
            }
        }

        // Write if there are faceFields
        if (faceFields.size() && Pstream::master())
        {
            os  << "CELL_DATA "
                << nFaces << nl
                << "FIELD attributes " << faceFields.size() << nl;
        }
        forAll(faceFields, fieldi)
        {
            const word& fieldName = faceFields[fieldi];
            FOR_ALL_FIELD_TYPES(WriteGeoField, VolField);
            FOR_ALL_FIELD_TYPES(WriteGeoField, SurfaceField);
        }

        // Write if there are pointFields
        if (pointFields.size() && Pstream::master())
        {
            os  << "POINT_DATA "
                << nPoints << nl
                << "FIELD attributes " << pointFields.size() << nl;
        }
        forAll(pointFields, fieldi)
        {
            const word& fieldName = pointFields[fieldi];
            FOR_ALL_FIELD_TYPES(WritePointField);
        }

        // Write the time series file
        if (Pstream::master())
        {
            timeSeries.writeTimeSeries(patchName);
        }
    }
}
