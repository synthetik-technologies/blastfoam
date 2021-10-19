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

Description
    GGI interpolation functions

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::interpolate
(
    const Field<Type>& ff,
    Field<Type>& result,
    const labelListList& addr,
    const scalarListList& weights
)
{
    forAll (result, faceI)
    {
        const labelList& curAddr = addr[faceI];
        const scalarList& curWeights = weights[faceI];

        result[faceI] = pTraits<Type>::zero;

        forAll (curAddr, i)
        {
            result[faceI] += ff[curAddr[i]]*curWeights[i];
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
(
    const Field<Type>& ff,
    Field<Type>& result,
    const labelList& mask,
    const labelListList& addr,
    const scalarListList& weights
)
{
    forAll (mask, maskI)
    {
        // Pick the masked face
        const label faceI = mask[maskI];

        const labelList& curAddr = addr[faceI];

        const scalarList& curWeights = weights[faceI];

        // Clear condensed list entry: masked faces only
        result[maskI] = pTraits<Type>::zero;

        forAll (curAddr, i)
        {
            // Put the result into condensed list: masked faces only
            result[maskI] += ff[curAddr[i]]*curWeights[i];
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::setFaces
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff,
    const labelList& facesToSet
)
{
    // Loop throuh fully uncovered faces and set the coefficients
    forAll (facesToSet, ftsI)
    {
        // Get face index
        const label& faceI = facesToSet[ftsI];

        // Set field for this face
        ff[faceI] = fieldToSet[faceI];
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedSetFaces
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff,
    const labelList& mask,
    const labelList& facesToSet
)
{
    // Note: tricky algorithm
    // In order for a face to be set it needs to be both in the
    // mask and in selection of faces that are set (addr).
    // This implies an n-squared search, but we can avoid it by
    // using the fact that both lists are ordered.

    label maskAddrI = 0;

    forAll (facesToSet, ftsI)
    {
        // Pick the face index
        const label& faceI = facesToSet[ftsI];

        // Search through the mask
        for (; maskAddrI < mask.size(); ++maskAddrI)
        {
            if (faceI == mask[maskAddrI])
            {
                // Found masked face, set the field
                ff[maskAddrI] = fieldToSet[maskAddrI];

                break;
            }
            else if (mask[maskAddrI] > faceI)
            {
                // Gone beyond my index: my face is not present in the mask
                // Go one back and check for next uncovered face
                if (maskAddrI > 0)
                {
                    --maskAddrI;
                }

                break;
            }
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::scalePartial
(
    Field<Type>& ff,
    const labelList& partiallyUncoveredAddr,
    const scalarField& uncoveredFractions
)
{
    // Loop through partially covered faces and scale them up
    forAll (partiallyUncoveredAddr, pcfI)
    {
        ff[partiallyUncoveredAddr[pcfI]] /=
            (1.0 - uncoveredFractions[pcfI]);
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedScalePartial
(
    Field<Type>& ff,
    const labelList& mask,
    const labelList& partiallyUncoveredAddr,
    const scalarField& uncoveredFractions
)
{
    // Note: tricky algorithm
    // In order for a face to be set it needs to be both in the
    // mask and in selection of faces that are set (addr).
    // This implies an n-squared search, but we can avoid it by
    // using the fact that both lists are ordered.

    label maskAddrI = 0;

    forAll (partiallyUncoveredAddr, pcfI)
    {
        // Pick partially covered face
        const label& faceI = partiallyUncoveredAddr[pcfI];

        for (; maskAddrI < mask.size(); ++maskAddrI)
        {
            if (faceI == mask[maskAddrI])
            {
                // Found masked partially covered face, scale it up
                ff[maskAddrI] /= (1.0 - uncoveredFractions[pcfI]);

                break;
            }
            else if (mask[maskAddrI] > faceI)
            {
                // Gone beyond my index: my face is not present in the mask
                // Go one back and check for next uncovered face
                if (maskAddrI > 0)
                {
                    --maskAddrI;
                }

                break;
            }
        }
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::addToPartialFaces
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff,
    const labelList& partiallyUncoveredAddr,
    const scalarField& uncoveredFractions
)
{
    // Loop through partially covered faces and add the weighted field
    forAll (partiallyUncoveredAddr, pcfI)
    {
        // Get face index
        const label& faceI = partiallyUncoveredAddr[pcfI];

        // Add to partially covered face
        ff[faceI] += uncoveredFractions[pcfI]*fieldToAdd[faceI];
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedAddToPartialFaces
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff,
    const labelList& mask,
    const labelList& partiallyUncoveredAddr,
    const scalarField& uncoveredFractions
)
{
    // Note: tricky algorithm
    // In order for a face to be partially covered it needs to be both in the
    // mask and in selection of faces that are partially covered. This implies
    // an n-squared search, but we can avoid it by using the fact that both
    // lists are ordered.

    label maskAddrI = 0;

    forAll (partiallyUncoveredAddr, pcfI)
    {
        // Pick partially covered face
        const label& faceI = partiallyUncoveredAddr[pcfI];

        for (; maskAddrI < mask.size(); ++maskAddrI)
        {
            if (faceI == mask[maskAddrI])
            {
                // Found masked partially covered face, add to it
                ff[maskAddrI] += uncoveredFractions[pcfI]*fieldToAdd[maskAddrI];

                break;
            }
            else if (mask[maskAddrI] > faceI)
            {
                // Gone beyond my index: my face is not present in the mask
                // Go one back and check for next uncovered face
                if (maskAddrI > 0)
                {
                    --maskAddrI;
                }

                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
GGIInterpolation<MasterPatch, SlavePatch>::masterToSlave
(
    const Field<Type>& ff
) const
{
    if (ff.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "GGIInterpolation::masterToSlave(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            slavePatch_.size(),
            pTraits<Type>::zero
        )
    );

    // Do interpolation
    Field<Type>& result = tresult.ref();

    if (this->doTransform() && pTraits<Type>::rank > 0)
    {
        // Transform master data to slave
        Field<Type> transformFF;

        if (reverseT_.size() == 1)
        {
            // Constant transform
            transformFF = transform(reverseT_[0], ff);
        }
        else
        {
            // Full patch transform
            transformFF = transform(reverseT_, ff);
        }

        GGIInterpolation<MasterPatch, SlavePatch>::interpolate
        (
            transformFF,
            result,
            this->slaveAddr(),
            this->slaveWeights()
        );
    }
    else
    {
        GGIInterpolation<MasterPatch, SlavePatch>::interpolate
        (
            ff,
            result,
            this->slaveAddr(),
            this->slaveWeights()
        );
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
GGIInterpolation<MasterPatch, SlavePatch>::masterToSlave
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = masterToSlave(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedMasterToSlave
(
    const Field<Type>& ff,
    Field<Type>& result,
    const labelList& mask
) const
{
    if (ff.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation::maskedMasterToSlave\n"
            "(\n"
            "    const Field<Type>& ff,\n"
            "    Field<Type>& result,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    if (result.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation::maskedMasterToSlave\n"
            "(\n"
            "    const Field<Type>& ff,\n"
            "    Field<Type>& result,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "result field does not correspond to mask. Field size: "
            << result.size() << " mask size: " << mask.size()
            << abort(FatalError);
    }

    if (this->doTransform() && pTraits<Type>::rank > 0)
    {
        // Transform master data to slave
        Field<Type> transformFF(ff.size());

        if (reverseT_.size() == 1)
        {
            // Constant transform
            // Transform only masked elements.  HJ, 25/May/2016
            forAll (mask, maskI)
            {
                transformFF[mask[maskI]] =
                    transform(reverseT_[0], ff[mask[maskI]]);
            }
        }
        else
        {
            // Full patch transform
            // Transform only masked elements.  HJ, 25/May/2016
            forAll (mask, maskI)
            {
                transformFF[mask[maskI]] =
                    transform(reverseT_[mask[maskI]], ff[mask[maskI]]);
            }
        }

        GGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
        (
            transformFF,
            result,
            mask,
            this->slaveAddr(),
            this->slaveWeights()
        );
    }
    else
    {
        GGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
        (
            ff,
            result,
            mask,
            this->slaveAddr(),
            this->slaveWeights()
        );
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
GGIInterpolation<MasterPatch, SlavePatch>::slaveToMaster
(
    const Field<Type>& ff
) const
{
    if (ff.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "GGIInterpolation::slaveToMaster"
            "(const Field<Type> ff)"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            masterPatch_.size(),
            pTraits<Type>::zero
        )
    );

    // Do interpolation
    Field<Type>& result = tresult.ref();

    if (this->doTransform() && pTraits<Type>::rank > 0)
    {
        // Transform slave data to master
        Field<Type> transformFF;
        if (forwardT_.size() == 1)
        {
            // Constant transform
            transformFF = transform(forwardT_[0], ff);
        }
        else
        {
            // Full patch transform
            transformFF = transform(forwardT_, ff);
        }

        GGIInterpolation<MasterPatch, SlavePatch>::interpolate
        (
            transformFF,
            result,
            this->masterAddr(),
            this->masterWeights()
        );
    }
    else
    {
        GGIInterpolation<MasterPatch, SlavePatch>::interpolate
        (
            ff,
            result,
            this->masterAddr(),
            this->masterWeights()
        );
    }

    return tresult;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
tmp<Field<Type> >
GGIInterpolation<MasterPatch, SlavePatch>::slaveToMaster
(
    const tmp<Field<Type> >& tff
) const
{
    tmp<Field<Type> > tint = slaveToMaster(tff());
    tff.clear();
    return tint;
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedSlaveToMaster
(
    const Field<Type>& ff,
    Field<Type>& result,
    const labelList& mask
) const
{
    if (ff.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation::maskedSlaveToMaster"
            "(\n"
            "    const Field<Type>& ff,\n"
            "    Field<Type>& result,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size() << " field size: " << ff.size()
            << abort(FatalError);
    }

    if (result.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation::maskedSlaveToMaster\n"
            "(\n"
            "    const Field<Type>& ff,\n"
            "    Field<Type>& result,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "result field does not correspond to mask. Field size: "
            << result.size() << " mask size: " << mask.size()
            << abort(FatalError);
    }

    if (this->doTransform() && pTraits<Type>::rank > 0)
    {
        // Transform slave data to master
        Field<Type> transformFF(ff.size());
        if (forwardT_.size() == 1)
        {
            // Constant transform
            // Transform only masked elements.  HJ, 25/May/2016
            forAll (mask, maskI)
            {
                transformFF[mask[maskI]] =
                    transform(forwardT_[0], ff[mask[maskI]]);
            }
        }
        else
        {
            // Full patch transform
            // Transform only masked elements.  HJ, 25/May/2016
            forAll (mask, maskI)
            {
                transformFF[mask[maskI]] =
                    transform(forwardT_[mask[maskI]], ff[mask[maskI]]);
            }
        }

        GGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
        (
            transformFF,
            result,
            mask,
            this->masterAddr(),
            this->masterWeights()
        );
    }
    else
    {
        GGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
        (
            ff,
            result,
            mask,
            this->masterAddr(),
            this->masterWeights()
        );
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::setUncoveredFacesMaster
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff
) const
{
    if
    (
        (fieldToSet.size() != masterPatch_.size())
     || (ff.size() != masterPatch_.size())
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "setUncoveredFacesMaster\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << abort(FatalError);
    }

    setFaces
    (
        fieldToSet,
        ff,
        uncoveredMasterFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedSetUncoveredFacesMaster
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedSetUncoveredFacesMaster\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << masterPatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedSetFaces
    (
        fieldToSet,
        ff,
        mask,
        uncoveredMasterFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::setUncoveredFacesSlave
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff
) const
{
    if
    (
        (fieldToSet.size() != slavePatch_.size())
     || (ff.size() != slavePatch_.size())
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "setUncoveredFacesSlave\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << abort(FatalError);
    }

    setFaces
    (
        fieldToSet,
        ff,
        uncoveredSlaveFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedSetUncoveredFacesSlave
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedSetUncoveredFacesSlave\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff\n,"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << slavePatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedSetFaces
    (
        fieldToSet,
        ff,
        mask,
        uncoveredSlaveFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::setPartialFacesMaster
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff
) const
{
    if
    (
        (fieldToSet.size() != masterPatch_.size())
     || (ff.size() != masterPatch_.size())
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "setPartialFacesMaster\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << abort(FatalError);
    }

    setFaces
    (
        fieldToSet,
        ff,
        partiallyUncoveredMasterFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedSetPartialFacesMaster
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedSetPartialFacesMaster\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << masterPatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedSetFaces
    (
        fieldToSet,
        ff,
        mask,
        partiallyUncoveredMasterFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::setPartialFacesSlave
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff
) const
{
    if
    (
        (fieldToSet.size() != slavePatch_.size())
     || (ff.size() != slavePatch_.size())
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "setPartialFacesSlave\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << abort(FatalError);
    }

    setFaces
    (
        fieldToSet,
        ff,
        partiallyUncoveredSlaveFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedSetPartialFacesSlave
(
    const Field<Type>& fieldToSet,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedSetPartialFacesSlave\n"
            "(\n"
            "    const Field<Type>& fieldToSet\n,"
            "    Field<Type>& ff\n,"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << slavePatch_.size()
            << " field to set size: " << fieldToSet.size()
            << " field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedSetFaces
    (
        fieldToSet,
        ff,
        mask,
        partiallyUncoveredSlaveFaces()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::scalePartialMaster
(
    Field<Type>& ff
) const
{
    if (ff.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "scalePartialMaster\n"
            "(\n"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size()
            << " scale field size: " << ff.size()
            << abort(FatalError);
    }

    scalePartial
    (
        ff,
        partiallyUncoveredMasterFaces(),
        masterFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedScalePartialMaster
(
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedScalePartialMaster\n"
            "(\n"
            "    Field<Type>& ff,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << masterPatch_.size()
            << " scale field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedScalePartial
    (
        ff,
        mask,
        partiallyUncoveredMasterFaces(),
        masterFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::scalePartialSlave
(
    Field<Type>& ff
) const
{
    if (ff.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "scalePartialSlave\n"
            "(\n"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size()
            << " scale field size: " << ff.size()
            << abort(FatalError);
    }

    scalePartial
    (
        ff,
        partiallyUncoveredSlaveFaces(),
        slaveFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedScalePartialSlave
(
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedScalePartialSlave\n"
            "(\n"
            "    Field<Type>& ff\n,"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << slavePatch_.size()
            << " scale field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedScalePartial
    (
        ff,
        mask,
        partiallyUncoveredSlaveFaces(),
        slaveFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::addToPartialFacesMaster
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff
) const
{
    if
    (
        fieldToAdd.size() != masterPatch_.size()
     || ff.size() != masterPatch_.size()
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "addToPartialFacesMaster\n"
            "(\n"
            "    const Field<Type>& fieldToAdd\n,"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size()
            << " field to add size: " << fieldToAdd.size()
            << " scale field size: " << ff.size()
            << abort(FatalError);
    }

    addToPartialFaces
    (
        fieldToAdd,
        ff,
        partiallyUncoveredMasterFaces(),
        masterFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedAddToPartialFacesMaster
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedAddToPartialFacesMaster\n"
            "(\n"
            "    const Field<Type>& fieldToAdd\n,"
            "    Field<Type>& ff,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << masterPatch_.size()
            << " field to set size: " << fieldToAdd.size()
            << " scale field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedAddToPartialFaces
    (
        fieldToAdd,
        ff,
        mask,
        partiallyUncoveredMasterFaces(),
        masterFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::addToPartialFacesSlave
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff
) const
{
    if
    (
        fieldToAdd.size() != slavePatch_.size()
     || ff.size() != slavePatch_.size()
    )
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "addToPartialFacesSlave\n"
            "(\n"
            "    const Field<Type>& fieldToAdd\n,"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size()
            << " field to add size: " << fieldToAdd.size()
            << " scale field size: " << ff.size()
            << abort(FatalError);
    }

    addToPartialFaces
    (
        fieldToAdd,
        ff,
        partiallyUncoveredSlaveFaces(),
        slaveFaceUncoveredFractions()
    );
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void GGIInterpolation<MasterPatch, SlavePatch>::maskedAddToPartialFacesSlave
(
    const Field<Type>& fieldToAdd,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void GGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedAddToPartialFacesSlave\n"
            "(\n"
            "    const Field<Type>& fieldToAdd\n,"
            "    Field<Type>& ff\n,"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch (mask) size: "
            << slavePatch_.size()
            << " field to add size: " << fieldToAdd.size()
            << " scale field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedAddToPartialFaces
    (
        fieldToAdd,
        ff,
        mask,
        partiallyUncoveredSlaveFaces(),
        slaveFaceUncoveredFractions()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
