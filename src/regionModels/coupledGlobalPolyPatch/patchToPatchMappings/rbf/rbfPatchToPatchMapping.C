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

void rbfPatchToPatchMapping::makeZoneAToZoneBInterpolator() const
{
    if (zoneAToZoneBInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!"
            << abort(FatalError);
    }

    Info<< "Create RBF interpolator from " << globalPatchA().patchName()
        << " to " << globalPatchB().patchName() << endl;

    const vectorField& zoneAFaceCentres = zoneA().faceCentres();
    const vectorField& zoneBFaceCentres = zoneB().faceCentres();

    zoneAToZoneBInterpolatorPtr_ =
        autoPtr<RBFInterpolation>
        (
            new RBFInterpolation
            (
                RBFFunctions::TPS::typeName,
                dict_,
                zoneAFaceCentres,
                zoneBFaceCentres
            )
        );

    // Check interpolation error
    vectorField zoneAFaceCentresAtZoneB
    (
        zoneAToZoneBInterpolatorPtr_->interpolate(zoneAFaceCentres)
    );
    const scalar maxDist = gMax
    (
        mag(zoneAFaceCentresAtZoneB - zoneBFaceCentres)
    );

    Info<< "    face interpolation error: " << maxDist << endl;
}


const RBFInterpolation&
rbfPatchToPatchMapping::zoneAToZoneBInterpolator() const
{
    if (!zoneAToZoneBInterpolatorPtr_.valid())
    {
        makeZoneAToZoneBInterpolator();
    }

    return zoneAToZoneBInterpolatorPtr_();
}


void rbfPatchToPatchMapping::makeZoneBToZoneAInterpolator() const
{
    if (zoneBToZoneAInterpolatorPtr_.valid())
    {
        FatalErrorInFunction
            << "Pointer already set!"
            << abort(FatalError);
    }

    Info<< "Create RBF interpolator from " << globalPatchB().patchName()
        << " to " << globalPatchA().patchName() << endl;

    const vectorField& zoneAFaceCentres = zoneA().faceCentres();
    const vectorField& zoneBFaceCentres = zoneB().faceCentres();

    zoneBToZoneAInterpolatorPtr_ =
        autoPtr<RBFInterpolation>
        (
            new RBFInterpolation
            (
                RBFFunctions::TPS::typeName,
                dict_,
                zoneBFaceCentres,
                zoneAFaceCentres
            )
        );

    // Check interpolation error
    vectorField zoneBFaceCentresAtZoneA
    (
        zoneBToZoneAInterpolatorPtr_->interpolate(zoneBFaceCentres)
    );
    const scalar maxDist = gMax
    (
        mag(zoneBFaceCentresAtZoneA - zoneAFaceCentres)
    );

    Info<< "    face interpolation error: " << maxDist << endl;
}


const RBFInterpolation&
rbfPatchToPatchMapping::zoneBToZoneAInterpolator() const
{
    if (!zoneBToZoneAInterpolatorPtr_.valid())
    {
        makeZoneBToZoneAInterpolator();
    }

    return zoneBToZoneAInterpolatorPtr_();
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
    zoneAToZoneBInterpolatorPtr_(NULL),
    zoneBToZoneAInterpolatorPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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
