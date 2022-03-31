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
#include "scalarMatrices.H"
namespace Foam
{
namespace patchToPatchMappings
{

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void rbfPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    Info<< "Interpolating face values using RBF" << endl;

    // Check field sizes are correct
    patchToPatchMapping::checkFieldSizes
    (
        fromZone.size(), toZone.size(), fromField.size(), toField.size()
    );

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is zoneA; toZone is zoneB
        zoneAToZoneBInterpolator().interpolate(fromField, toField);
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is zoneA; fromZone is zoneB
        zoneBToZoneAInterpolator().interpolate(fromField, toField);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown zones!" << abort(FatalError);
    }
}


template<class Type>
void rbfPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    Info<< "Interpolating point values using RBF" << endl;

    // Check field sizes are correct
    patchToPatchMapping::checkFieldSizes
    (
        fromZone.nPoints(), toZone.nPoints(), fromField.size(), toField.size()
    );


    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is zoneA; toZone is zoneB
        zoneAToZoneBInterpolator().interpolate(fromField, toField);
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is zoneA; fromZone is zoneB
        zoneBToZoneAInterpolator().interpolate(fromField, toField);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown zones!" << abort(FatalError);
    }
}

// ************************************************************************* //

} // End namespace patchToPatchMappings

} // End namespace Foam


// ************************************************************************* //
