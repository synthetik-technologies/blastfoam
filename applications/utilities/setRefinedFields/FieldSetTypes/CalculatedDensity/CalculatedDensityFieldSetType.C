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

#include "CalculatedDensityFieldSetType.H"
#include "calcAngleFraction.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetTypes::CalculatedDensity<Type, Patch, Mesh>::CalculatedDensity
(
    const fvMesh& mesh,
    const word& fieldName,
    const labelList& selectedCells,
    Istream& is,
    const bool write
)
:
    FieldSetType<Type, Patch, Mesh>(mesh, fieldName, selectedCells, is, write),
    volInt_(pTraits<Type>(is))
{
    if (this->good_)
    {
        scalar V = 0.0;
        forAll(selectedCells, i)
        {
            V += mesh.V()[selectedCells[i]];
        }
        reduce(V, sumOp<scalar>());
        V /= calcAngleFraction(mesh);
        if (V > small)
        {
            this->good_ = true;
            density_ = volInt_/V;
            setField();
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::FieldSetTypes::CalculatedDensity<Type, Patch, Mesh>::~CalculatedDensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
void Foam::FieldSetTypes::CalculatedDensity<Type, Patch, Mesh>::setField()
{
    if (this->selectedCells_.size() == this->mesh_.nCells())
    {
        (*this->fieldPtr_).primitiveFieldRef() = density_;
    }
    else
    {
        forAll(this->selectedCells_, i)
        {
            (*this->fieldPtr_)[this->selectedCells_[i]] = density_;
        }
    }
    FieldSetType<Type, Patch, Mesh>::setField();
}

