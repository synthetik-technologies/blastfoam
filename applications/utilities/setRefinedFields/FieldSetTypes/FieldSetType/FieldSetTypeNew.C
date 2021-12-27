/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

template<class Type, template<class> class Patch, class Mesh>
Foam::autoPtr<Foam::FieldSetType<Type, Patch, Mesh>>
Foam::FieldSetType<Type, Patch, Mesh>::New
(
    const word& fieldDesc,
    const word& fieldName,
    const fvMesh& mesh,
    const dictionary& dict,
    const labelList& selectedCells,
    Istream& is,
    const bool write
)
{
    word fieldSetType(fieldDesc);
    word fieldType(FieldType::typeName);

    fieldSetType.replaceAll(fieldType, word::null);

    typename dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTablePtr_->find(fieldSetType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown fieldSetType type " << fieldSetType
            << " for fieldSetType " << FieldType::typeName
            << "Valid fieldSetType types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    return cstrIter()(mesh, dict, fieldName, selectedCells, is, write);
}
