/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Store field maxes over all times
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

#include "fieldMax.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::fieldMax::createMax
(
    const word& fieldName,
    const word& maxFieldName
)
{
    if (!obr_.foundObject<Type>(fieldName))
    {
        return;
    }

    Log << "    Reading/initialising field " << maxFieldName << endl;

    if (obr_.foundObject<Type>(maxFieldName))
    {}
    else if (obr_.found(maxFieldName))
    {
        Log << "    Cannot allocate average field " << maxFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;
    }
    else if (obr_.found(fieldName))
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);

        // Store on registry
        obr_.store
        (
            new Type
            (
                IOobject
                (
                    maxFieldName,
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


template<class Type>
void Foam::functionObjects::fieldMax::createOld
(
    const word& maxFieldName,
    PtrListDictionary<Type>& old
)
{
    if (obr_.foundObject<Type>(maxFieldName))
    {
        const Type& f = obr_.lookupObject<Type>(maxFieldName);

        // Store unregistered fields so fields are not updated with refinement
        old.resize(old.size() + 1);
        old.set
        (
            old.size()-1,
            maxFieldName,
            new Type
            (
                IOobject
                (
                    maxFieldName,
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


template<class Type>
void Foam::functionObjects::fieldMax::updateMax
(
    const word& fieldName,
    const word& maxFieldName,
    PtrListDictionary<Type>& old
)
{
    if (obr_.foundObject<Type>(fieldName))
    {
        const Type& baseField = obr_.lookupObject<Type>(fieldName);
        Type& field = obr_.lookupObjectRef<Type>(maxFieldName);
        field = max(field, baseField);

        if (cellMap_)
        {
            const labelList& cellMap = *cellMap_;
            const labelList& rCellMap = *rCellMap_;

            const Type& fOld = old[maxFieldName];

            forAll(cellMap, i)
            {
                label celli = cellMap[i];
                if (celli > -1)
                {
                    field[i] = max(field[i], fOld[celli]);
                }
            }

            forAll(rCellMap, i)
            {
                label index = rCellMap[i];

                if (index < -1)
                {
                    label celli = -index-2;
                    field[celli] = max(field[celli], fOld[i]);
                }
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldMax::mapMax
(
    const word& fieldName,
    const mapPolyMesh& meshMap,
    const PtrListDictionary<Type>& old
)
{
    if (obr_.foundObject<Type>(fieldName))
    {
        Type& f = obr_.lookupObjectRef<Type>(fieldName);

        if (cellMap_)
        {
            const labelList& cellMap = *cellMap_;
            const labelList& rCellMap = *rCellMap_;

            const Type& fOld = old[fieldName];
            forAll(cellMap, i)
            {
                label celli = cellMap[i];
                if (celli > -1)
                {
                    f[i] = max(f[i], fOld[celli]);
                }
            }

            forAll(rCellMap, i)
            {
                label index = rCellMap[i];

                if (index < -1)
                {
                    label celli = -index-2;
                    f[celli] = max(f[celli], fOld[i]);
                }
            }
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldMax::writeField(const word& fieldName)
{
    if (obr_.foundObject<Type>(fieldName))
    {
        obr_.lookupObject<Type>(fieldName).write();
    }
}

// ************************************************************************* //
