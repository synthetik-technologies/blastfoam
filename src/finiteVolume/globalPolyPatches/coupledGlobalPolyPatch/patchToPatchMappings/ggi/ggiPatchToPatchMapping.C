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

#include "ggiPatchToPatchMapping.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace patchToPatchMappings
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ggiPatchToPatchMapping, 0);
addToRunTimeSelectionTable
(
    patchToPatchMapping, ggiPatchToPatchMapping, dictionary
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void ggiPatchToPatchMapping::makeInterpolator() const
{
    if (interpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer is already set!"
            << abort(FatalError);
    }

    if (debug)
    {
        Info<< "Create GGI zone-to-zone interpolator" << endl;
    }

    interpolatorPtr_.set
    (
        new GGIInterpolation<standAlonePatch, standAlonePatch>
        (
            zoneA(),
            zoneB(),
            tensorField(0),
            tensorField(0),
            vectorField(0), // Slave-to-master separation. Bug fix
            true,           // Patch data is complete on all processors
            SMALL,          // Non-overlapping face tolerances
            SMALL,
            true,           // Rescale weighting factors
            ggiInterpolation::BB_OCTREE
        )
    );

    if (debug)
    {
        checkZoneAToZoneBError();
        checkZoneBToZoneAError();

        Info<< "Number of uncovered master faces: "
            << interpolatorPtr_().uncoveredMasterFaces().size() << nl
            << "Number of uncovered slave faces: "
            << interpolatorPtr_().uncoveredSlaveFaces().size() << nl << endl;
    }
}

const GGIInterpolation<standAlonePatch, standAlonePatch>&
ggiPatchToPatchMapping::interpolator() const
{
    if (interpolatorPtr_.empty())
    {
        makeInterpolator();
    }

    return interpolatorPtr_();
}


void ggiPatchToPatchMapping::checkZoneAToZoneBError() const
{
    // Reference to patch face centres
    const vectorField& patchAFaceCentres = patchA().faceCentres();

    // Construct global zone field
    const vectorField zoneAFaceCentres
    (
        globalPatchA().patchFaceToGlobal(patchAFaceCentres)
    );

    // Interpolate global zone field from A to B
    const vectorField zoneBFaceCentres
    (
        interpolator().masterToSlave(zoneAFaceCentres)
    );

    // Extract local patch field
    const vectorField patchBFaceCentres
    (
        globalPatchB().globalFaceToPatch(zoneBFaceCentres)
    );

    // Print maximum error
    Info<< "interface-to-interface face error: "
        << gMax(mag(patchBFaceCentres - patchB().faceCentres()))
        << endl;
}


void ggiPatchToPatchMapping::checkZoneBToZoneAError() const
{
    const vectorField zoneBPointsAtFluid
    (
        interpolator().slaveToMasterPointInterpolate
        (
            zoneB().localPoints()
        )
    );

    const vectorField& zoneAPoints = zoneA().localPoints();

    Info<< "interface-to-interface point error: "
        << gMax(mag(zoneAPoints - zoneBPointsAtFluid))
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

ggiPatchToPatchMapping::ggiPatchToPatchMapping
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
    interpolatorPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void ggiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferFaces<scalar>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferFaces<vector>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferFaces<symmTensor>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferFaces<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<tensor>& fromField,  // from field
    Field<tensor>& toField           // to field
) const
{
    transferFaces<tensor>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<scalar>& fromField,  // from field
    Field<scalar>& toField           // to field
) const
{
    transferPoints<scalar>(fromZone, toZone, fromField, toField);
}

void ggiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<vector>& fromField,  // from field
    Field<vector>& toField           // to field
) const
{
    transferPoints<vector>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<symmTensor>& fromField,  // from field
    Field<symmTensor>& toField           // to field
) const
{
    transferPoints<symmTensor>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<sphericalTensor>& fromField,  // from field
    Field<sphericalTensor>& toField           // to field
) const
{
    transferPoints<sphericalTensor>(fromZone, toZone, fromField, toField);
}


void ggiPatchToPatchMapping::transferPoints
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
