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

#include "amiPatchToPatchMapping.H"

namespace Foam
{
namespace patchToPatchMappings
{

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //

template<class Type>
void amiPatchToPatchMapping::transferFaces
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    if (debug)
    {
        Info<< "Interpolating face values using AMI" << endl;
    }

    // Check field sizes are correct
    patchToPatchMapping::checkFieldSizes
    (
        fromZone.size(), toZone.size(), fromField.size(), toField.size()
    );

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is the master (zoneA); toZone is the slave (zoneB)
        toField = interpolator().interpolateToTarget(fromField);
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is the master (zoneA); fromZone is the slave (zoneB)
        toField = interpolator().interpolateToSource(fromField);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown zones!" << abort(FatalError);
    }
}


template<class Type>
void amiPatchToPatchMapping::transferPoints
(
    const standAlonePatch& fromZone, // from zone
    const standAlonePatch& toZone,   // to zone
    const Field<Type>& fromField,    // from field
    Field<Type>& toField             // to field
) const
{
    if (debug)
    {
        Info<< "Interpolating point values using AMI" << endl;
    }

    // Check field sizes are correct
    patchToPatchMapping::checkFieldSizes
    (
        fromZone.nPoints(), toZone.nPoints(), fromField.size(), toField.size()
    );

    // Check if fromZone is zoneA or zoneB by checking the memory address
    if (&fromZone == &zoneA() && &toZone == &zoneB())
    {
        // fromZone is the master (zoneA); toZone is the slave (zoneB)
        //toField = interpolator().interpolateToTargetPoints(fromField);

        // Escape the interpolation if there are no faces in the target patch
        if (zoneB().nPoints() == 0)
        {
            return;
        }

        const List<face>& zoneAFaces = zoneA().localFaces();
        const List<labelPair>& addr = zoneBPointAddr();
        const FieldField<Field, double>& weights = zoneBPointWeights();

        forAll(toField, pointI)
        {
            if (addr[pointI].first() > -1)
            {
                const face& hitFace = zoneAFaces[addr[pointI].first()];
                const label pI = addr[pointI].second();
                const Type ctrF = average(Field<Type>(fromField, hitFace));

                toField[pointI] =
                    weights[pointI][0]*fromField[hitFace[pI]]
                  + weights[pointI][1]*fromField[hitFace.nextLabel(pI)]
                  + weights[pointI][2]*ctrF;
            }
            else
            {
                FatalErrorInFunction
                    << "zoneB point addressing is not correct"
                    << abort(FatalError);
            }
        }
    }
    else if (&toZone == &zoneA() && &fromZone == &zoneB())
    {
        // toZone is the master (zoneA); fromZone is the slave (zoneB)
        //toField = interpolateToSourcePoints(interpolator(), fromField);

        // Escape the interpolation if there are no faces in the target patch
        if (zoneA().nPoints() == 0)
        {
            return;
        }

        const List<face>& zoneBFaces = zoneB().localFaces();
        const List<labelPair>& addr = zoneAPointAddr();
        const FieldField<Field, double>& weights = zoneAPointWeights();

        forAll(toField, pointI)
        {
            if (addr[pointI].first() > -1)
            {
                const face& hitFace = zoneBFaces[addr[pointI].first()];

                const label pI = addr[pointI].second();

                const Type ctrF = average(Field<Type>(fromField, hitFace));

                toField[pointI] =
                    weights[pointI][0]*fromField[hitFace[pI]]
                  + weights[pointI][1]*fromField[hitFace.nextLabel(pI)]
                  + weights[pointI][2]*ctrF;
            }
            else
            {
                FatalErrorInFunction
                    << "zoneA point addressing is not correct"
                    << abort(FatalError);
            }
        }
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
