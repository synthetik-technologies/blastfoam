/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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

#include "energyImmersedBoundaryFvPatchScalarField.H"
#include "immersedBoundaryFvPatchFields.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BaseType>
Foam::energyImmersedBoundaryFvPatchScalarField<BaseType>::
energyImmersedBoundaryFvPatchScalarField
(
    volScalarField& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
:
    BaseType(f, dict, ibo, "T"),
    mesh_(f.mesh())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BaseType>
void Foam::energyImmersedBoundaryFvPatchScalarField<BaseType>::updateCoeffs() const
{
    BaseType::updateCoeffs();
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName
            (
                basicThermo::dictName,
                this->field_.group()
            )
        );

    const labelList patchCells(this->ibm_.patchExternalCells());
    labelList cells(patchCells.size());
    Field<scalar> values(this->values_);
    labelList map(patchCells.size(), -1);

    label I = 0;
    forAll(patchCells, i)
    {
        label celli = patchCells[i];
        if (celli >= 0)
        {
            map[I] = i;
            values[I] = values[i];
            cells[I] = celli;
            I++;
        }
    }
    if (I == 0)
    {
        return;
    }
    map.resize(I);
    cells.resize(I);
    values.resize(I);
    values = thermo.he(values, cells);
    this->values_ = Zero;
    forAll(map, i)
    {
        this->values_[map[i]] = values[i];
    }
}


template<class BaseType>
void Foam::energyImmersedBoundaryFvPatchScalarField<BaseType>::setValues()
{
    BaseType::updateCoeffs();
//     const basicThermo& thermo =
//         mesh_.lookupObject<basicThermo>
//         (
//             IOobject::groupName
//             (
//                 basicThermo::dictName,
//                 this->field_.group()
//             )
//         );
//     scalar avg(sum(this->values_*this->ibm_.magSf())/sum(this->ibm_.magSf()));
//
//     const labelList& faceCells(this->ibm_.patchExternalCells());
//     this->ibm_.setInternal
//     (
//         this->field_,
//         thermo.he(Field<scalar>(faceCells.size(), avg), faceCells)()
//     );
}

// ************************************************************************* //
