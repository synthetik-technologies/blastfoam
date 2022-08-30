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

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::CalculatedDensity<Type, FSType>::CalculatedDensity
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
    volInt_(pTraits<Type>(is))
{
    if (this->good_)
    {
        scalar V = 0.0;
        forAll(selectedIndices, i)
        {
            V += mesh.V()[selectedIndices[i]];
        }
        reduce(V, sumOp<scalar>());
        V /= calcAngleFraction(mesh);
        if (V > small)
        {
            this->good_ = true;
            density_ = volInt_/V;
            this->setField();
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
Foam::FieldSetTypes::CalculatedDensity<Type, FSType>::~CalculatedDensity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class FSType>
void
Foam::FieldSetTypes::CalculatedDensity<Type, FSType>::getInternalField
(
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    f = density_;
}


template<class Type, template<class> class FSType>
void
Foam::FieldSetTypes::CalculatedDensity<Type, FSType>::getBoundaryField
(
    const label patchi,
    const labelList& indices,
    const UIndirectList<vector>& pts,
    UIndirectList<Type>& f
)
{
    f = density_;
}
