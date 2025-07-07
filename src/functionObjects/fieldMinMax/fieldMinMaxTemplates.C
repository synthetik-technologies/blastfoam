/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies:  |  Store field maxes over all times
07-09-2021 Synthetik Applied Technologies:  |  Added minimum option
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "fieldMinMax.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
Type Foam::functionObjects::fieldMinMax::compute
(
    const Type& a,
    const Type& b
) const
{
    Type modA(a);
    Type modB(b);
    switch (mode_)
    {
        case modeType::mag:
        {
            scalar sqrtNCmpts =
                sqrt(scalar(pTraits<Type>::nComponents));

            modA = mag(a)/sqrtNCmpts*pTraits<Type>::one;
            modA = mag(b)/sqrtNCmpts*pTraits<Type>::one;
            break;
        }
        case modeType::cmptMag:
        {
            modA = cmptMag(a);
            modA = cmptMag(b);
            break;
        }
        default:
            break;
    }
    return
        minMax_ == minMaxType::min
      ? min(modA, modB)
      : max(modA, modB);
}


template<class FieldType>
void Foam::functionObjects::fieldMinMax::compute
(
    FieldType& res,
    const FieldType& a,
    const FieldType& b
) const
{
    forAll(res, i)
    {
        res[i] = compute(a[i], b[i]);
    }
}


template<class FieldType>
bool Foam::functionObjects::fieldMinMax::createMinMax
(
    const word& fieldName
)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        return true;
    }
    if (!obr_.foundObject<FieldType>(fieldName))
    {
        return false;
    }

    typedef typename FieldType::value_type Type;

    Info << "    Reading/initialising field " << computeFieldName << endl;
    const FieldType& baseField = obr_.lookupObject<FieldType>(fieldName);

    // Store on registry
    FieldType* fPtr =
        new FieldType
        (
            IOobject
            (
                computeFieldName,
                obr_.time().timeName(obr_.time().startTime().value()),
                obr_,
                restartOnRestart_
              ? IOobject::NO_READ
              : IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            baseField.mesh(),
            dimensioned<Type>(baseField.dimensions(), Zero)
        );
    if (!fPtr->headerOk())
    {
        (*fPtr) = baseField;
    }

    fPtr->store(fPtr);

    return true;
}


template<class FieldType>
bool Foam::functionObjects::fieldMinMax::storeOld
(
    const word& fieldName
)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        const FieldType& f = obr_.lookupObject<FieldType>(computeFieldName);

        if (!oldFields_.found(computeFieldName))
        {
            // Store unregistered fields so fields are not updated with refinement
            oldFields_.set
            (
                computeFieldName,
                new FieldType
                (
                    IOobject
                    (
                        computeFieldName,
                        obr_.time().name(),
                        obr_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    f
                )
            );
        }
        else
        {
            dynamicCast<FieldType&>(*oldFields_[computeFieldName]).reset(f);
        }
        return true;
    }
    return false;
}


template<class FieldType>
bool Foam::functionObjects::fieldMinMax::update
(
    const word& fieldName
)
{
    if (obr_.foundObject<FieldType>(fieldName))
    {
        const word computeFieldName(computedName(fieldName));
        objectRegistry::const_iterator oldIter =
            obr_.find(computeFieldName);
        if (oldIter == obr_.cend())
        {
            createMinMax<FieldType>(fieldName);
        }

        const FieldType& baseField = obr_.lookupObject<FieldType>(fieldName);
        FieldType& field = obr_.lookupObjectRef<FieldType>(computeFieldName);
        compute
        (
            field.primitiveFieldRef(),
            field.primitiveField(),
            baseField.primitiveField()
        );
        forAll(field.boundaryField(), patchi)
        {
            compute
            (
                field.boundaryFieldRef()[patchi],
                field.boundaryField()[patchi],
                baseField.boundaryField()[patchi]
            );
        }
        return true;
    }
    return false;
}


template<class FieldType>
bool Foam::functionObjects::fieldMinMax::mapField
(
    const word& fieldName,
    const polyTopoChangeMap& map
)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        FieldType& f = obr_.lookupObjectRef<FieldType>(computeFieldName);

        const labelList& cellMap = map.cellMap();
        const labelList& rCellMap = map.reverseCellMap();

        HashPtrTable<regIOobject>::iterator fOldIter =
            oldFields_.find(computeFieldName);
        const FieldType& fOld =
            *dynamic_cast<const FieldType*>(fOldIter());
        forAll(cellMap, i)
        {
            label celli = cellMap[i];
            if (celli > -1)
            {
                f[i] = compute(f[i], fOld[celli]);
            }
        }

        forAll(rCellMap, i)
        {
            label index = rCellMap[i];

            if (index < -1)
            {
                label celli = -index-2;
                f[celli] = compute(f[celli], fOld[i]);
            }
        }

        return true;
    }
    return false;
}


// ************************************************************************* //
