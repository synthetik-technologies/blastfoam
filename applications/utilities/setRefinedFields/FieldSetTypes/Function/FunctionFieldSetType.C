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

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetTypes::Function<Type, Patch, Mesh>::Function
(
    const fvMesh& mesh,
    const word& fieldName,
    const labelList& selectedCells,
    Istream& is,
    const bool write
)
:
    FieldSetType<Type, Patch, Mesh>(mesh, fieldName, selectedCells, is, write),
    func_()
{
    dictionary dict;
    dict.add(fieldName, dictionary(is));

    func_ = Function3<Type>::New(fieldName, dict);

    if (this->good_)
    {
        setField();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetTypes::Function<Type, Patch, Mesh>::~Function()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::FieldSetTypes::Function<Type, Patch, Mesh>::setField()
{
    Info<<(*this->fieldPtr_).name()<<endl;
    const vectorField& C = this->mesh_.C();
    forAll(this->selectedCells_, i)
    {
        label celli = this->selectedCells_[i];
        (*this->fieldPtr_)[celli] =
            func_->value
            (
                C[celli].x(),
                C[celli].y(),
                C[celli].z()
            );
    }

    FieldSetType<Type, Patch, Mesh>::setField();
}

