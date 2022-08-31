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

#include "FunctionFieldSetType.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::Function<Type, FSType>::Function
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& fieldName,
    const labelList& selectedIndices,
    Istream& is,
    const bool write
)
:
    FSType<Type>
    (
        mesh,
        dict,
        fieldName,
        selectedIndices,
        is,
        write
    ),
    func_()
{
    dictionary funcDict(dict, dictionary());
    funcDict.add(fieldName, dictionary(is));

    func_ = Function3<Type>::New(fieldName, funcDict);

    if (this->good_)
    {
        this->setField();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::Function<Type, FSType>::~Function()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
void Foam::FieldSetTypes::Function<Type, FSType>::getInternalField
(
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    forAll(indices, i)
    {
        f[i] =
            func_->value
            (
                pts[i].x(),
                pts[i].y(),
                pts[i].z()
            );
    }
}


template<class Type, template<class> class FSType>
void Foam::FieldSetTypes::Function<Type, FSType>::getBoundaryField
(
    const label patchi,
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    forAll(indices, i)
    {
        f[i] =
            func_->value
            (
                pts[i].x(),
                pts[i].y(),
                pts[i].z()
            );
    }
}
