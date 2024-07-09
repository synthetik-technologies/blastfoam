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

#include "patchToPatchMapping.H"

// * * * * * * * * * * * * * Static Members Functions  * * * * * * * * * * * //

template<class SubListA, class SubListB>
inline void Foam::patchToPatchMapping::transferListList
(
    List<SubListA>& a,
    List<SubListB>& b
)
{
    a.setSize(b.size());
    forAll(a, i)
    {
        a[i].transfer(b[i]);
    }
}


template<class Type>
inline void Foam::patchToPatchMapping::rDistributeListList
(
    const label size,
    const mapDistribute& map,
    List<List<Type>>& data
)
{
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>(),
        size,
        map.constructMap(),
        false,
        map.subMap(),
        false,
        data,
        ListAppendEqOp<Type>(),
        flipOp(),
        List<Type>()
    );
}


template<class Type>
inline void Foam::patchToPatchMapping::rDistributeListList
(
    const label size,
    const mapDistribute& map,
    List<DynamicList<Type>>& data
)
{
    List<List<Type>> tdata;
    transferListList(tdata, data);
    rDistributeListList(size, map, tdata);
    transferListList(data, tdata);
}


template<class Type, class LabelList, class ScalarList>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::interpolate
(
    const List<LabelList>& localOtherData,
    const List<ScalarList>& weights,
    const autoPtr<mapDistribute>& otherMapPtr,
    const Field<Type>& otherFld
)
{
    tmp<Field<Type>> tlocalOtherFld;
    if (otherMapPtr.valid())
    {
        tlocalOtherFld = tmp<Field<Type>>(new Field<Type>(otherFld));
        otherMapPtr->distribute(tlocalOtherFld.ref());
    }
    const Field<Type>& localOtherFld =
        tlocalOtherFld.valid() ? tlocalOtherFld() : otherFld;

    tmp<Field<Type>> tfld
    (
        new Field<Type>
        (
            localOtherData.size(),
            pTraits<Type>::one*Foam::NaN
        )
    );
    Field<Type>& fld = tfld.ref();

    forAll(localOtherData, datai)
    {
        const labelList& otherData = localOtherData[datai];
        if (otherData.size())
        {
            scalar sumW = 0;
            Type sumWData = Zero;

            const ScalarList& ws = weights[datai];
            forAll(otherData, i)
            {
                const scalar w = ws[i];
                sumW += w;
                sumWData += localOtherFld[otherData[i]]*w;
            }
            fld[datai] = sumWData/sumW;
        }
    }
    return tfld;
}


template<class Type, class LabelList, class ScalarList>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::interpolate
(
    const List<LabelList>& localOtherData,
    const List<ScalarList>& weights,
    const autoPtr<mapDistribute>& otherMapPtr,
    const Field<Type>& otherFld,
    const Field<Type>& leftOverFld
)
{
    tmp<Field<Type>> tlocalOtherFld;
    if (otherMapPtr.valid())
    {
        tlocalOtherFld = tmp<Field<Type>>(new Field<Type>(otherFld));
        otherMapPtr->distribute(tlocalOtherFld.ref());
    }
    const Field<Type>& localOtherFld =
        tlocalOtherFld.valid() ? tlocalOtherFld() : otherFld;

    tmp<Field<Type>> tfld
    (
        new Field<Type>
        (
            localOtherData.size(),
            pTraits<Type>::one*Foam::NaN
        )
    );
    Field<Type>& fld = tfld.ref();

    forAll(localOtherData, datai)
    {
        const labelList& otherData = localOtherData[datai];
        if (otherData.size())
        {
            scalar sumW = 0;
            Type sumWData(Zero);

            const ScalarList& ws = weights[datai];
            forAll(otherData, i)
            {
                const scalar w = ws[i];
                sumW += w;
                sumWData += localOtherFld[otherData[i]]*w;
            }
            fld[datai] = sumWData + (1.0 - sumW)*leftOverFld[datai];
        }
    }
    return tfld;
}

// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::transferToTgt
(
    const Field<Type>& srcField,
    const Field<Type>& unmapped
) const
{
    if (!isNull(unmapped))
    {
        return interpolate
        (
            localSrcToTgt_,
            tgtWeights()(),
            srcMapPtr_,
            srcField,
            unmapped
        );
    }
    else
    {
        return interpolate
        (
            localSrcToTgt_,
            tgtWeights()(),
            srcMapPtr_,
            srcField
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::patchToPatchMapping::transferToTgt
(
    const tmp<Field<Type>>& tsrcField,
    const Field<Type>& unmapped
) const
{
    return transferToTgt(tsrcField(), unmapped);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>  Foam::patchToPatchMapping::transferToSrc
(
    const Field<Type>& tgtField,
    const Field<Type>& unmapped
) const
{
    if (!isNull(unmapped))
    {
        return interpolate
        (
            localTgtToSrc_,
            srcWeights()(),
            tgtMapPtr_,
            tgtField,
            unmapped
        );
    }
    else
    {
        return interpolate
        (
            localTgtToSrc_,
            srcWeights()(),
            tgtMapPtr_,
            tgtField
        );
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>  Foam::patchToPatchMapping::transferToSrc
(
    const tmp<Field<Type>>& ttgtField,
    const Field<Type>& unmapped
) const
{
    return transferToSrc(ttgtField(), unmapped);
}


// ************************************************************************* //
