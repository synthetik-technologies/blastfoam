/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "ggiPolyPatch.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::fastExpand
(
    const UList<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    // with communication
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::fastExpand\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )   << "Incorrect patch field size.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (localParallel())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::fastExpand"
            "("
            "    const Field<Type>& ff"
            ") const"
        )   << "Requested expand on local parallel.  This is not allowed"
            << abort(FatalError);
    }

    // Replaced old comms algorithm.  HJ, 31/May/2016
    // HJ, 4/Jun/2011

    // HR, 10/Jul/2013
    // This function requires send-receive-addressing, but usage is not
    // symmetric across processors. Hence trigger re-calculate at this point
    if (Pstream::parRun() && !localParallel())
    {
        map();
        shadow().map();
    }

    // New version: mapDistribute

    if (Pstream::parRun())
    {
        // Optimised mapDistribute

        // Prepare return field: expand the field to zone size
        tmp<Field<Type> > texpandField
        (
            new Field<Type>(ff)
        );
        Field<Type>& expandField = texpandField.ref();

        map().distribute(expandField);

        return texpandField;
    }
    else
    {
        // Serial.  Expand the field to zone size

        tmp<Field<Type> > texpandField
        (
            new Field<Type>(zone().size())  // filled with nans
        );
        Field<Type>& expandField = texpandField.ref();

        const labelList& zAddr = zoneAddressing();

        forAll (zAddr, i)
        {
            expandField[zAddr[i]] = ff[i];
        }

        return texpandField;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::interpolate
(
    const Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != shadow().size())
    {
        FatalErrorIn
        (
            "tmp<Field<Type> > ggiPolyPatch::interpolate\n"
            "(\n"
            "    const Field<Type>& ff\n"
            ") const"
        )   << "Incorrect slave patch field size.  Field size: "
            << ff.size() << " patch size: " << shadow().size()
            << abort(FatalError);
    }

    // New.  HJ, 12/Jun/2011
    if (localParallel())
    {
        // No expansion or filtering needed.  HJ, 4/Jun/2011

        if (empty())
        {
            // Patch empty, no interpolation
            return tmp<Field<Type> >(new Field<Type>());
        }

        // Interpolate field
        if (owner())
        {
            return patchToPatch().slaveToMaster(ff);
        }
        else
        {
            return patchToPatch().masterToSlave(ff);
        }
    }
    else
    {
        // Expand shadow
        Field<Type> expandField(shadow().fastExpand(ff));

        tmp<Field<Type> > tresult(new Field<Type>(size()));
        Field<Type>& result = tresult.ref();

        if (owner())
        {
            patchToPatch().maskedSlaveToMaster
            (
                expandField,
                result,
                zoneAddressing()
            );
        }
        else
        {
            patchToPatch().maskedMasterToSlave
            (
                expandField,
                result,
                zoneAddressing()
            );
        }

        return tresult;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::ggiPolyPatch::interpolate
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = interpolate(tff());
    tff.clear();
    return tint;
}


template<class Type>
void Foam::ggiPolyPatch::setUncoveredFaces
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "template<class Type> void ggiPolyPatch::setUncoveredFaces\n"
            "(\n"
            "    const Field<Type>& fieldToSet,\n"
            "    Field<Type>& ff\n"
            ") const"
        )   << "Incorrect patch field size for setting.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (bridgeOverlap())
    {
        if (empty())
        {
            // Patch empty, no bridging
            return;
        }

        if (localParallel())
        {
            if (owner())
            {
                patchToPatch().setUncoveredFacesMaster(fieldToSet, ff);
            }
            else
            {
                patchToPatch().setUncoveredFacesSlave(fieldToSet, ff);
            }
        }
        else
        {
            // Note: since bridging is only a local operation
            if (owner())
            {
                patchToPatch().maskedSetUncoveredFacesMaster
                (
                    fieldToSet,
                    ff,
                    zoneAddressing()
                );
            }
            else
            {
                patchToPatch().maskedSetUncoveredFacesSlave
                (
                    fieldToSet,
                    ff,
                    zoneAddressing()
                );
            }
        }
    }
}


template<class Type>
void Foam::ggiPolyPatch::setPartialFaces
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != size())
    {
        FatalErrorInFunction
            << "Incorrect patch field size for setting.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (bridgeOverlap())
    {
        if (empty())
        {
            // Patch empty, no bridging
            return;
        }

        if (localParallel())
        {
            if (owner())
            {
                patchToPatch().setPartialFacesMaster(fieldToSet, ff);
            }
            else
            {
                patchToPatch().setPartialFacesSlave(fieldToSet, ff);
            }
        }
        else
        {
            // Note: since bridging is only a local operation
            if (owner())
            {
                patchToPatch().maskedSetPartialFacesMaster
                (
                    fieldToSet,
                    ff,
                    zoneAddressing()
                );
            }
            else
            {
                patchToPatch().maskedSetPartialFacesSlave
                (
                    fieldToSet,
                    ff,
                    zoneAddressing()
                );
            }
        }
    }
}


template<class Type>
void Foam::ggiPolyPatch::scalePartialFaces(Field<Type>& ff) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "template<class Type> ggiPolyPatch::scalePartialFaces\n"
            "(\n"
            "    Field<Type>& ff,\n"
            ") const"
        )   << "Incorrect patch field size for scaling.  Field size: "
            << ff.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (bridgeOverlap())
    {
        if (empty())
        {
            // Patch empty, no bridging
            return;
        }

        if (localParallel())
        {
            if (owner())
            {
                patchToPatch().scalePartialMaster(ff);
            }
            else
            {
                patchToPatch().scalePartialSlave(ff);
            }
        }
        else
        {
            // Note: since bridging is only a local operation
            if (owner())
            {
                patchToPatch().maskedScalePartialMaster
                (
                    ff,
                    zoneAddressing()
                );
            }
            else
            {
                patchToPatch().maskedScalePartialSlave
                (
                    ff,
                    zoneAddressing()
                );
            }
        }
    }
}


template<class Type>
void Foam::ggiPolyPatch::addToPartialFaces
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff
) const
{
    // Check and expand the field from patch size to zone size
    if (ff.size() != size())
    {
        FatalErrorIn
        (
            "template<class Type> ggiPolyPatch::addToPartialFaces\n"
            "(\n"
            "    const Field<Type>& fieldToAdd,\n"
            "    Field<Type>& ff,\n"
            ") const"
        )   << "Incorrect patch field size for adding.  Field size: "
            << ff.size() << " field to add size: "
            << fieldToAdd.size() << " patch size: " << size()
            << abort(FatalError);
    }

    if (bridgeOverlap())
    {
        if (empty())
        {
            // Patch empty, no bridging
            return;
        }

        if (localParallel())
        {
            if (owner())
            {
                patchToPatch().addToPartialFacesMaster(fieldToAdd, ff);
            }
            else
            {
                patchToPatch().addToPartialFacesSlave(fieldToAdd, ff);
            }
        }
        else
        {
            // Note: since bridging is only a local operation
            if (owner())
            {
                patchToPatch().maskedAddToPartialFacesMaster
                (
                    fieldToAdd,
                    ff,
                    zoneAddressing()
                );
            }
            else
            {
                patchToPatch().maskedAddToPartialFacesSlave
                (
                    fieldToAdd,
                    ff,
                    zoneAddressing()
                );
            }
        }
    }
}


// ************************************************************************* //
