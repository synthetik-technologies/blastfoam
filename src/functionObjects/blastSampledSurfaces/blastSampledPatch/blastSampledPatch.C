/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "blastSampledPatch.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blastSampledSurfaces
{
    defineTypeNameAndDebug(patch, 0);
    addToRunTimeSelectionTable(blastSampledSurface, patch, word);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blastSampledSurfaces::patch::patch
(
    const word& name,
    const polyMesh& mesh,
    const wordReList& patchNames,
    const bool triangulate
)
:
    blastSampledSurface(name, mesh),
    patchNames_(patchNames),
    triangulate_(triangulate),
    needsUpdate_(true)
{}


Foam::blastSampledSurfaces::patch::patch
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    blastSampledSurface(name, mesh, dict),
    patchNames_(dict.lookup("patches")),
    triangulate_(dict.lookupOrDefault("triangulate", false)),
    needsUpdate_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blastSampledSurfaces::patch::~patch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::blastSampledSurfaces::patch::patchIDs() const
{
    if (patchIDs_.empty())
    {
        patchIDs_ = mesh().boundaryMesh().patchSet
        (
            patchNames_,
            false
        ).sortedToc();
    }
    return patchIDs_;
}


bool Foam::blastSampledSurfaces::patch::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::blastSampledSurfaces::patch::expire()
{
    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    blastSampledSurface::clearGeom();
    MeshStorage::clear();
    patchIDs_.clear();
    patchIndex_.clear();
    patchFaceLabels_.clear();
    patchStart_.clear();

    needsUpdate_ = true;
    return true;
}


bool Foam::blastSampledSurfaces::patch::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    label sz = 0;
    label np = 0;
    forAll(patchIDs(), i)
    {
        label patchi = patchIDs()[i];
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        if (isA<emptyPolyPatch>(pp))
        {
            FatalErrorInFunction
                << "Cannot sample an empty patch. Patch " << pp.name()
                << exit(FatalError);
        }

        sz += pp.size();
        np += pp.nPoints();
    }

    // For every face (or triangle) the originating patch and local face in the
    // patch.
    patchIndex_.setSize(sz);
    patchFaceLabels_.setSize(sz);
    patchStart_.setSize(patchIDs().size());
    labelList meshFaceLabels(sz);

    patchPointIndex_.setSize(np);
    patchPointLabels_.setSize(np);
    patchPointStart_.setSize(patchIDs().size());

    sz = 0;
    np = 0;
    forAll(patchIDs(), i)
    {
        label patchi = patchIDs()[i];

        patchStart_[i] = sz;
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        forAll(pp, j)
        {
            patchIndex_[sz] = i;
            patchFaceLabels_[sz] = j;
            meshFaceLabels[sz] = pp.start()+j;
            sz++;
        }

        patchPointStart_[i] = np;
        const labelList& meshPoints = pp.meshPoints();
        forAll(meshPoints, j)
        {
            patchPointIndex_[np] = i;
            patchPointLabels_[np] = j;
            np++;
        }
    }

    indirectPrimitivePatch allPatches
    (
        IndirectList<face>(mesh().faces(), meshFaceLabels),
        mesh().points()
    );

    this->storedPoints() = allPatches.localPoints();
    this->storedFaces()  = allPatches.localFaces();


    // triangulate uses remapFaces()
    // - this is somewhat less efficient since it recopies the faces
    // that we just created, but we probably don't want to do this
    // too often anyhow.
    if (triangulate_)
    {
        MeshStorage::triangulate();
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    needsUpdate_ = false;
    return true;
}


// remap action on triangulation
void Foam::blastSampledSurfaces::patch::remapFaces(const labelUList& faceMap)
{
    // recalculate the cells cut
    if (notNull(faceMap) && faceMap.size())
    {
        MeshStorage::remapFaces(faceMap);
        patchFaceLabels_ = labelList
        (
            UIndirectList<label>(patchFaceLabels_, faceMap)
        );
        patchIndex_ = labelList
        (
            UIndirectList<label>(patchIndex_, faceMap)
        );

        // Redo patchStart.
        if (patchIndex_.size() > 0)
        {
            patchStart_[patchIndex_[0]] = 0;
            for (label i = 1; i < patchIndex_.size(); i++)
            {
                if (patchIndex_[i] != patchIndex_[i-1])
                {
                    patchStart_[patchIndex_[i]] = i;
                }
            }
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::blastSampledSurfaces::patch::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField> Foam::blastSampledSurfaces::patch::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::blastSampledSurfaces::patch::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField> Foam::blastSampledSurfaces::patch::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField> Foam::blastSampledSurfaces::patch::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField> Foam::blastSampledSurfaces::patch::sample
(
    const surfaceScalarField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::vectorField> Foam::blastSampledSurfaces::patch::sample
(
    const surfaceVectorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::blastSampledSurfaces::patch::sample
(
    const surfaceSphericalTensorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::symmTensorField> Foam::blastSampledSurfaces::patch::sample
(
    const surfaceSymmTensorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::tensorField> Foam::blastSampledSurfaces::patch::sample
(
    const surfaceTensorField& sField
) const
{
    return sampleField(sField);
}


Foam::tmp<Foam::scalarField> Foam::blastSampledSurfaces::patch::sample
(
    const pointScalarField& pField
) const
{
    return sampleField(pField);
}


Foam::tmp<Foam::vectorField> Foam::blastSampledSurfaces::patch::sample
(
    const pointVectorField& pField
) const
{
    return sampleField(pField);
}


Foam::tmp<Foam::sphericalTensorField> Foam::blastSampledSurfaces::patch::sample
(
    const pointSphericalTensorField& pField
) const
{
    return sampleField(pField);
}


Foam::tmp<Foam::symmTensorField> Foam::blastSampledSurfaces::patch::sample
(
    const pointSymmTensorField& pField
) const
{
    return sampleField(pField);
}


Foam::tmp<Foam::tensorField> Foam::blastSampledSurfaces::patch::sample
(
    const pointTensorField& pField
) const
{
    return sampleField(pField);
}


Foam::tmp<Foam::scalarField> Foam::blastSampledSurfaces::patch::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::blastSampledSurfaces::patch::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::sphericalTensorField> Foam::blastSampledSurfaces::patch::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::blastSampledSurfaces::patch::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::blastSampledSurfaces::patch::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::blastSampledSurfaces::patch::print(Ostream& os) const
{
    os  << "patch: " << name() << " :"
        << "  patches:" << patchNames()
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
