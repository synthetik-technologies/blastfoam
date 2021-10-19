/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
void newGGIInterpolation<MasterPatch, SlavePatch>::interpolate
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
void newGGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
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
void newGGIInterpolation<MasterPatch, SlavePatch>::bridge
(
    const Field<Type>& bridgeField,
    Field<Type>& ff,
    const labelList& addr
)
{
    forAll (addr, faceI)
    {
        ff[addr[faceI]] = bridgeField[addr[faceI]];
    }
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void newGGIInterpolation<MasterPatch, SlavePatch>::maskedBridge
(
    const Field<Type>& bridgeField,
    Field<Type>& ff,
    const labelList& mask,
    const labelList& uncoveredFaces
)
{
    // Note: tricky algorithm
    // In order for a face to be bridged it needs to be both in the
    // mask and in selection of faces that are bridged (addr).
    // This implies an n-squared search, but we can avoid it by
    // using the fact that both lists are ordered.

    label maskAddrI = 0;

    forAll (uncoveredFaces, uncoI)
    {
        // Pick the uncovered face
        const label faceI = uncoveredFaces[uncoI];

        // Search through the mask
        for (; maskAddrI < mask.size(); maskAddrI++)
        {
            if (faceI == mask[maskAddrI])
            {
                // Found masked bridged face
                // Put the result into condensed list: masked faces only
                ff[maskAddrI] = bridgeField[maskAddrI];

                break;
            }
            else if (mask[maskAddrI] > faceI)
            {
                // Gone beyond my index: my face is not present in the mask
                // Go one back and check for next uncovered face
                if (maskAddrI > 0)
                {
                    maskAddrI--;
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
newGGIInterpolation<MasterPatch, SlavePatch>::masterToSlave
(
    const Field<Type>& ff
) const
{
    if (ff.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "newGGIInterpolation::masterToSlave(const Field<Type> ff)"
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
    Field<Type>& result = tresult();

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

        newGGIInterpolation<MasterPatch, SlavePatch>::interpolate
        (
            transformFF,
            result,
            this->slaveAddr(),
            this->slaveWeights()
        );
    }
    else
    {
        newGGIInterpolation<MasterPatch, SlavePatch>::interpolate
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
newGGIInterpolation<MasterPatch, SlavePatch>::masterToSlave
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
void newGGIInterpolation<MasterPatch, SlavePatch>::maskedMasterToSlave
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
            "void newGGIInterpolation::maskedMasterToSlave\n"
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
            "void newGGIInterpolation::maskedMasterToSlave\n"
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

        newGGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
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
        newGGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
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
newGGIInterpolation<MasterPatch, SlavePatch>::slaveToMaster
(
    const Field<Type>& ff
) const
{
    if (ff.size() != slavePatch_.size())
    {
        FatalErrorIn
        (
            "newGGIInterpolation::slaveToMaster"
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
    Field<Type>& result = tresult();

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

        newGGIInterpolation<MasterPatch, SlavePatch>::interpolate
        (
            transformFF,
            result,
            this->masterAddr(),
            this->masterWeights()
        );
    }
    else
    {
        newGGIInterpolation<MasterPatch, SlavePatch>::interpolate
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
newGGIInterpolation<MasterPatch, SlavePatch>::slaveToMaster
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
void newGGIInterpolation<MasterPatch, SlavePatch>::maskedSlaveToMaster
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
            "void newGGIInterpolation::maskedSlaveToMaster"
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
            "void newGGIInterpolation::maskedSlaveToMaster\n"
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

        newGGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
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
        newGGIInterpolation<MasterPatch, SlavePatch>::maskedInterpolate
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
void newGGIInterpolation<MasterPatch, SlavePatch>::bridgeMaster
(
    const Field<Type>& bridgeField,
    Field<Type>& ff
) const
{
    if
    (
        bridgeField.size() != masterPatch_.size()
     || ff.size() != masterPatch_.size())
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::bridgeMaster\n"
            "(\n"
            "    const Field<Type>& bridgeField,\n"
            "    Field<Type>& ff\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size()
            << " bridge field size: " << bridgeField.size()
            << " field size: " << ff.size()
            << abort(FatalError);
    }

    bridge(bridgeField, ff, uncoveredMasterFaces());
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void newGGIInterpolation<MasterPatch, SlavePatch>::maskedBridgeMaster
(
    const Field<Type>& bridgeField,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedBridgeMaster\n"
            "(\n"
            "    const Field<Type>& bridgeField,\n"
            "    Field<Type>& ff,\n"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << masterPatch_.size()
            << " bridge field size: " << bridgeField.size()
            << " field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedBridge(bridgeField, ff, mask, uncoveredMasterFaces());
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void newGGIInterpolation<MasterPatch, SlavePatch>::bridgeSlave
(
    const Field<Type>& bridgeField,
    Field<Type>& ff
) const
{
    if
    (
        bridgeField.size() != slavePatch_.size()
     || ff.size() != slavePatch_.size()
    )
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "bridgeSlave\n"
            "(\n"
            "    const Field<Type>& bridgeField,\n"
            "    Field<Type>& ff"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size()
            << " bridge field size: " << bridgeField.size()
            << " field size: " << ff.size()
            << abort(FatalError);
    }

    bridge(bridgeField, ff, uncoveredSlaveFaces());
}


template<class MasterPatch, class SlavePatch>
template<class Type>
void newGGIInterpolation<MasterPatch, SlavePatch>::maskedBridgeSlave
(
    const Field<Type>& bridgeField,
    Field<Type>& ff,
    const labelList& mask
) const
{
    if (ff.size() != mask.size())
    {
        FatalErrorIn
        (
            "void newGGIInterpolation<MasterPatch, SlavePatch>::"
            "maskedBridgeSlave\n"
            "(\n"
            "    const Field<Type>& bridgeField,\n"
            "    Field<Type>& ff\n,"
            "    const labelList& mask\n"
            ") const"
        )   << "given field does not correspond to patch. Patch size: "
            << slavePatch_.size()
            << " bridge field size: " << bridgeField.size()
            << " field size: " << ff.size()
            << " mask size: " << mask.size()
            << abort(FatalError);
    }

    maskedBridge(bridgeField, ff, mask, uncoveredSlaveFaces());
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
