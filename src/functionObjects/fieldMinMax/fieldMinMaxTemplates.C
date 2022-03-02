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
void Foam::functionObjects::fieldMinMax::createMinMax
(
    const word& fieldName
)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        return;
    }
    if (!obr_.foundObject<FieldType>(fieldName))
    {
        return;
    }

    Log << "    Reading/initialising field " << computeFieldName << endl;

    if (obr_.found(computeFieldName))
    {
        Log << "    Cannot allocate average field " << computeFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;
    }
    else if (obr_.found(fieldName))
    {
        const FieldType& baseField = obr_.lookupObject<FieldType>(fieldName);

        // Store on registry
        obr_.store
        (
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
                baseField,
                baseField.boundaryField()
            )
        );
    }
}


template<class FieldType>
void Foam::functionObjects::fieldMinMax::createOld
(
    const word& fieldName,
    HashPtrTable<FieldType>& old
)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        const FieldType& f = obr_.lookupObject<FieldType>(computeFieldName);

        // Store unregistered fields so fields are not updated with refinement
        old.insert
        (
            computeFieldName,
            new FieldType
            (
                IOobject
                (
                    computeFieldName,
                    obr_.time().timeName(),
                    obr_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                f
            )
        );
    }
}


template<class FieldType>
void Foam::functionObjects::fieldMinMax::update
(
    const word& fieldName,
    const HashPtrTable<FieldType>& old
)
{
    const word computeFieldName(computedName(fieldName));
    if (!old.found(computeFieldName))
    {
        createMinMax<FieldType>(fieldName);
    }

    if (obr_.foundObject<FieldType>(fieldName))
    {
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

        if (cellMap_.valid())
        {
            const labelList& cellMap = cellMap_();
            const labelList& rCellMap = rCellMap_();

            const FieldType& fOld = *old[computeFieldName];

            forAll(cellMap, i)
            {
                label celli = cellMap[i];
                if (celli > -1)
                {
                    field[i] = compute(field[i], fOld[celli]);
                }
            }

            forAll(rCellMap, i)
            {
                label index = rCellMap[i];

                if (index < -1)
                {
                    label celli = -index-2;
                    field[celli] = compute(field[celli], fOld[i]);
                }
            }
        }
    }
}


template<class FieldType>
void Foam::functionObjects::fieldMinMax::map
(
    const word& fieldName,
    const mapPolyMesh& meshMap,
    const HashPtrTable<FieldType>& old
)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        FieldType& f = obr_.lookupObjectRef<FieldType>(computeFieldName);

        if (cellMap_.valid())
        {
            const labelList& cellMap = cellMap_();
            const labelList& rCellMap = rCellMap_();

            const FieldType& fOld = *old[computeFieldName];
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
        }
    }
}


template<class FieldType>
void Foam::functionObjects::fieldMinMax::writeField(const word& fieldName)
{
    const word computeFieldName(computedName(fieldName));
    if (obr_.foundObject<FieldType>(computeFieldName))
    {
        writeObject(computeFieldName);
    }
}

// ************************************************************************* //
