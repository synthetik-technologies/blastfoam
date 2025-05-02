/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "blastSampledSurfaces.H"
#include "PatchTools.H"
#include "polyTopoChangeMap.H"
#include "OSspecific.H"
#include "writeFile.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(blastSampledSurfaces, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        blastSampledSurfaces,
        dictionary
    );
}
}

bool Foam::functionObjects::blastSampledSurfaces::verbose_ = false;
Foam::scalar Foam::functionObjects::blastSampledSurfaces::mergeTol_ = 1e-10;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::blastSampledSurfaces::writeGeometry() const
{
    // Write to time directory under outputPath_
    // Skip surface without faces (eg, a failed cut-plane)

    const fileName outputDir = outputPath_/mesh_.time().name();

    forAll(*this, surfI)
    {
        const blastSampledSurface& s = operator[](surfI);

        if (Pstream::parRun())
        {
            if (Pstream::master() && mergeList_[surfI].faces.size())
            {
                formatter_->write
                (
                    outputDir,
                    s.name(),
                    mergeList_[surfI].points,
                    mergeList_[surfI].faces
                );
            }
        }
        else if (s.faces().size())
        {
            formatter_->write
            (
                outputDir,
                s.name(),
                s.points(),
                s.faces()
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::blastSampledSurfaces::blastSampledSurfaces
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name, t),
    PtrList<blastSampledSurface>(),
    mesh_
    (
        refCast<const fvMesh>
        (
            t.lookupObject<objectRegistry>
            (
                dict.lookupOrDefault("region", polyMesh::defaultRegion)
            )
        )
    ),
    outputPath_(fileName::null),
    fieldSelection_(),
    interpolationScheme_(word::null),
    mergeList_(),
    formatter_(nullptr)
{
    outputPath_ =
        mesh_.time().globalPath()/functionObjects::writeFile::outputPrefix/name;

    if (mesh_.name() != fvMesh::defaultRegion)
    {
        outputPath_ = outputPath_/mesh_.name();
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::blastSampledSurfaces::~blastSampledSurfaces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::functionObjects::blastSampledSurfaces::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::functionObjects::blastSampledSurfaces::execute()
{
    return true;
}


bool Foam::functionObjects::blastSampledSurfaces::write()
{
    if (size())
    {
        // Finalise surfaces, merge points etc.
        update();

        // Create the output directory
        if (Pstream::master())
        {
            if (debug)
            {
                Pout<< "Creating directory "
                    << outputPath_/mesh_.time().name() << nl << endl;

            }

            mkDir(outputPath_/mesh_.time().name());
        }

        // Create a list of names of fields that are actually available
        wordList fieldNames;
        forAll(fields_, fieldi)
        {
            #define FoundFieldType(Type, nullArg)             \
              || foundObject<VolField<Type>>(fields_[fieldi]) \
              || foundObject<SurfaceField<Type>>(fields_[fieldi])
            if (false FOR_ALL_FIELD_TYPES(FoundFieldType))
            {
                fieldNames.append(fields_[fieldi]);
            }
            else
            {
                cannotFindObject(fields_[fieldi]);
            }
            #undef FoundFieldType
        }

        // Create table of cached interpolations, to prevent unnecessary work
        // when interpolating fields over multiple surfaces
        #define DeclareInterpolations(Type, nullArg) \
            HashPtrTable<interpolation<Type>> interpolation##Type##s;
        FOR_ALL_FIELD_TYPES(DeclareInterpolations);
        #undef DeclareInterpolations

        // Sample and write the surfaces
        forAll(*this, surfi)
        {
            const sampledSurface& s = operator[](surfi);

            #define GenerateFieldTypeValues(Type, nullArg) \
                PtrList<Field<Type>> field##Type##Values = \
                    sampleType<Type>(surfi, fieldNames, interpolation##Type##s);
            FOR_ALL_FIELD_TYPES(GenerateFieldTypeValues);
            #undef GenerateFieldTypeValues

            if (Pstream::parRun())
            {
                if
                (
                    Pstream::master()
                 && (mergeList_[surfi].faces.size() || writeEmpty_)
                )
                {
                    formatter_->write
                    (
                        outputPath_/mesh_.time().name(),
                        s.name(),
                        mergeList_[surfi].points,
                        mergeList_[surfi].faces,
                        fieldNames,
                        s.interpolate()
                        #define FieldTypeValuesParameter(Type, nullArg) \
                            , field##Type##Values
                        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                        #undef FieldTypeValuesParameter
                    );
                }
            }
            else
            {
                if (s.faces().size() || writeEmpty_)
                {
                    formatter_->write
                    (
                        outputPath_/mesh_.time().name(),
                        s.name(),
                        s.points(),
                        s.faces(),
                        fieldNames,
                        s.interpolate()
                        #define FieldTypeValuesParameter(Type, nullArg) \
                            , field##Type##Values
                        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
                        #undef FieldTypeValuesParameter
                    );
                }
            }
        }
    }

    return true;

    return true;
}


bool Foam::functionObjects::blastSampledSurfaces::read(const dictionary& dict)
{
    bool surfacesFound = dict.found("surfaces");

    if (surfacesFound)
    {
        dict.lookup("fields") >> fieldSelection_;

        dict.lookup("interpolationScheme") >> interpolationScheme_;
        const word writeType(dict.lookup("surfaceFormat"));

        // Define the surface formatter
        formatter_ = surfaceWriter::New(writeType, dict);

        PtrList<blastSampledSurface> newList
        (
            dict.lookup("surfaces"),
            blastSampledSurface::iNew(mesh_)
        );
        transfer(newList);

        if (Pstream::parRun())
        {
            mergeList_.setSize(size());
        }

        // Ensure all surfaces and merge information are expired
        expire();

        if (this->size())
        {
            Info<< "Reading surface description:" << nl;
            forAll(*this, surfI)
            {
                Info<< "    " << operator[](surfI).name() << nl;
            }
            Info<< endl;
        }
    }

    if (Pstream::master() && debug)
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample surfaces:" << nl << "(" << nl;

        forAll(*this, surfI)
        {
            Pout<< "  " << operator[](surfI) << endl;
        }
        Pout<< ")" << endl;
    }

    return true;
}


void Foam::functionObjects::blastSampledSurfaces::updateMesh
(
    const polyTopoChangeMap& mpm
)
{
    if (&mpm.mesh() == &mesh_)
    {
        expire();
    }

    // pointMesh and interpolation will have been reset in mesh.update
}


void Foam::functionObjects::blastSampledSurfaces::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        expire();
    }
}


void Foam::functionObjects::blastSampledSurfaces::readUpdate
(
    const polyMesh::readUpdateState state
)
{
    if (state != polyMesh::UNCHANGED)
    {
        expire();
    }
}


bool Foam::functionObjects::blastSampledSurfaces::needsUpdate() const
{
    forAll(*this, surfI)
    {
        if (operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::functionObjects::blastSampledSurfaces::expire()
{
    bool justExpired = false;

    forAll(*this, surfI)
    {
        if (operator[](surfI).expire())
        {
            justExpired = true;
        }

        // Clear merge information
        if (Pstream::parRun())
        {
            mergeList_[surfI].clear();
        }
    }

    // true if any surfaces just expired
    return justExpired;
}


bool Foam::functionObjects::blastSampledSurfaces::update()
{
    bool updated = false;

    if (!needsUpdate())
    {
        return updated;
    }

    // Serial: quick and easy, no merging required
    if (!Pstream::parRun())
    {
        forAll(*this, surfI)
        {
            if (operator[](surfI).update())
            {
                updated = true;
            }
        }

        return updated;
    }

    // Dimension as fraction of mesh bounding box
    scalar mergeDim = mergeTol_ * mesh_.bounds().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
            << mergeDim << " metre" << endl;
    }

    forAll(*this, surfI)
    {
        blastSampledSurface& s = operator[](surfI);

        if (s.update())
        {
            updated = true;
        }
        else
        {
            continue;
        }

        PatchTools::gatherAndMerge
        (
            mergeDim,
            primitivePatch
            (
                SubList<face>(s.faces(), s.faces().size()),
                s.points()
            ),
            mergeList_[surfI].points,
            mergeList_[surfI].faces,
            mergeList_[surfI].pointsMap
        );
    }

    return updated;
}


// ************************************************************************* //
