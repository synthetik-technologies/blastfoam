/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "rbfPatchToPatchMapping.H"
#include "addToRunTimeSelectionTable.H"
#include "TPSRBFFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rbfPatchToPatchMapping, 0);
addToRunTimeSelectionTable
(
    patchToPatchMapping,
    rbfPatchToPatchMapping,
    dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void rbfPatchToPatchMapping::makeZoneAToZoneBFaceInterpolator() const
{
    if (zoneAToZoneBFaceInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!"
            << abort(FatalError);
    }

    DebugInfo<< "Create RBF face interpolator from A to B, "
        << globalPatchA().patchName() << " to "
        << globalPatchB().patchName() << endl;

    const vectorField& zoneAFaceCentres = zoneA().faceCentres();
    const vectorField& zoneBFaceCentres = zoneB().faceCentres();

    zoneAToZoneBFaceInterpolatorPtr_ =
        autoPtr<RBFInterpolation>
        (
            new RBFInterpolation
            (
                dict_.lookupOrDefault("RBFFunction", RBFFunctions::TPS::typeName),
                dict_,
                zoneAFaceCentres,
                zoneBFaceCentres,
                CONSISTENT
            )
        );

    if (debug)
    {
        // Check interpolation error
        vectorField zoneAFaceCentresAtZoneB
        (
            zoneAToZoneBFaceInterpolatorPtr_->interpolate(zoneAFaceCentres)
        );
        const scalar maxDist = gMax
        (
            mag(zoneAFaceCentresAtZoneB - zoneBFaceCentres)
        );

        Info<< typeName << ": Face interpolation error= " << maxDist << endl;
    }
}


const RBFInterpolation&
rbfPatchToPatchMapping::zoneAToZoneBFaceInterpolator() const
{
    if (!zoneAToZoneBFaceInterpolatorPtr_.valid())
    {
        makeZoneAToZoneBFaceInterpolator();
    }

    return zoneAToZoneBFaceInterpolatorPtr_();
}


void rbfPatchToPatchMapping::makeZoneBToZoneAFaceInterpolator() const
{
    if (zoneBToZoneAFaceInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!"
            << abort(FatalError);
    }

    DebugInfo<< "Create RBF face interpolator from B to A, "
        << globalPatchB().patchName() << " to "
        << globalPatchA().patchName() << endl;

    const vectorField& zoneAFaceCentres = zoneA().faceCentres();
    const vectorField& zoneBFaceCentres = zoneB().faceCentres();

    zoneBToZoneAFaceInterpolatorPtr_ =
        autoPtr<RBFInterpolation>
        (
            new RBFInterpolation
            (
                dict_.lookupOrDefault("RBFFunction", RBFFunctions::TPS::typeName),
                dict_,
                zoneBFaceCentres,
                zoneAFaceCentres,
                CONSISTENT
            )
        );

    if (debug)
    {
        // Check interpolation error
        vectorField zoneBFaceCentresAtZoneA
        (
            zoneBToZoneAFaceInterpolatorPtr_->interpolate(zoneBFaceCentres)
        );
        const scalar maxDist = gMax
        (
            mag(zoneBFaceCentresAtZoneA - zoneAFaceCentres)
        );

        Info<< typeName << ": Face interpolation error= " << maxDist << endl;
    }
}


const RBFInterpolation&
rbfPatchToPatchMapping::zoneBToZoneAFaceInterpolator() const
{
    if (!zoneBToZoneAFaceInterpolatorPtr_.valid())
    {
        makeZoneBToZoneAFaceInterpolator();
    }

    return zoneBToZoneAFaceInterpolatorPtr_();
}


void rbfPatchToPatchMapping::makeZoneAToZoneBPointInterpolator() const
{
    if (zoneAToZoneBPointInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!"
            << abort(FatalError);
    }

    DebugInfo<< "Create RBF point interpolator from A to B, "
        << globalPatchA().patchName() << " to "
        << globalPatchB().patchName() << endl;

    const vectorField& zoneAPoints = zoneA().localPoints();
    const vectorField& zoneBPoints = zoneB().localPoints();

    zoneAToZoneBPointInterpolatorPtr_ =
        autoPtr<RBFInterpolation>
        (
            new RBFInterpolation
            (
                dict_.lookupOrDefault("RBFFunction", RBFFunctions::TPS::typeName),
                dict_,
                zoneAPoints,
                zoneBPoints,
                CONSISTENT
            )
        );

    if (debug)
    {
        // Check interpolation error
        vectorField zoneAPointsAtZoneB
        (
            zoneAToZoneBPointInterpolatorPtr_->interpolate(zoneAPoints)
        );
        const scalar maxDist = gMax
        (
            mag(zoneAPointsAtZoneB - zoneBPoints)
        );

        Info<< typeName << ": Point interpolation error= " << maxDist << endl;
    }
}


const RBFInterpolation&
rbfPatchToPatchMapping::zoneAToZoneBPointInterpolator() const
{
    if (!zoneAToZoneBPointInterpolatorPtr_.valid())
    {
        makeZoneAToZoneBPointInterpolator();
    }

    return zoneAToZoneBPointInterpolatorPtr_();
}


void rbfPatchToPatchMapping::makeZoneBToZoneAPointInterpolator() const
{
    if (zoneBToZoneAPointInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!"
            << abort(FatalError);
    }

    DebugInfo<< "Create RBF point interpolator from B to A, "
        << globalPatchB().patchName() << " to "
        << globalPatchA().patchName() << endl;

    const vectorField& zoneAPoints = zoneA().localPoints();
    const vectorField& zoneBPoints = zoneB().localPoints();

    zoneBToZoneAPointInterpolatorPtr_ =
        autoPtr<RBFInterpolation>
        (
            new RBFInterpolation
            (
                dict_.lookupOrDefault("RBFFunction", RBFFunctions::TPS::typeName),
                dict_,
                zoneBPoints,
                zoneAPoints,
                CONSISTENT
            )
        );

    if (debug)
    {
        // Check interpolation error
        vectorField zoneBPointsAtZoneA
        (
            zoneBToZoneAPointInterpolatorPtr_->interpolate(zoneBPoints)
        );
        const scalar maxDist = gMax
        (
            mag(zoneBPointsAtZoneA - zoneAPoints)
        );

        Info<< typeName << ": Point interpolation error= " << maxDist << endl;
    }
}


const RBFInterpolation&
rbfPatchToPatchMapping::zoneBToZoneAPointInterpolator() const
{
    if (!zoneBToZoneAPointInterpolatorPtr_.valid())
    {
        makeZoneBToZoneAPointInterpolator();
    }

    return zoneBToZoneAPointInterpolatorPtr_();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rbfPatchToPatchMapping::rbfPatchToPatchMapping
(
    const dictionary& dict,
    const primitivePatch& patchA,
    const primitivePatch& patchB,
    const globalPolyPatch& globalPatchA,
    const globalPolyPatch& globalPatchB
)
:
    patchToPatchMapping
    (
        typeName_(), dict, patchA, patchB, globalPatchA, globalPatchB
    ),
    dict_(dict),
    zoneAToZoneBFaceInterpolatorPtr_(NULL),
    zoneBToZoneAFaceInterpolatorPtr_(NULL),
    zoneAToZoneBPointInterpolatorPtr_(NULL),
    zoneBToZoneAPointInterpolatorPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

labelList rbfPatchToPatchMapping::unmappedFaces
(
    const standAlonePatch& patch
) const
{
    Field<scalar> res
    (
        (&patch == &(zoneA()))
      ? zoneBToZoneAFaceInterpolator().interpolate
        (
            scalarField(zoneB().size(), 1.0)
        )
      : zoneAToZoneBFaceInterpolator().interpolate
        (
            scalarField(zoneA().size(), 1.0)
        )
    );
    DynamicList<label> unmapped;
    forAll(res, i)
    {
        if (mag(res[i]) < small)
        {
            unmapped.append(i);
        }
    }
    return unmapped;
}

labelList rbfPatchToPatchMapping::unmappedPoints
(
    const standAlonePatch& patch
) const
{
    Field<scalar> res
    (
        (&patch == &(zoneA()))
      ? zoneBToZoneAPointInterpolator().interpolate
        (
            scalarField(zoneB().nPoints(), 1.0)
        )
      : zoneAToZoneBPointInterpolator().interpolate
        (
            scalarField(zoneA().nPoints(), 1.0)
        )
    );
    DynamicList<label> unmapped;
    forAll(res, i)
    {
        if (mag(res[i]) < small)
        {
            unmapped.append(i);
        }
    }
    return unmapped;
}

void rbfPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFaces<scalar>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFaces<vector>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferFaces<symmTensor>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferFaces<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferFaces<tensor>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPoints<scalar>(fromZone, toZone, fromField, toField);
}

void rbfPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferPoints<vector>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferPoints<symmTensor>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferPoints<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void rbfPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferPoints<tensor>(fromZone, toZone, fromField, toField);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchToPatchMappings

} // End namespace Foam

// ************************************************************************* //
