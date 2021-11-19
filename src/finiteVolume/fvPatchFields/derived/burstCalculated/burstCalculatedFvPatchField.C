/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "burstCyclicFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::burstCyclicFvPatchField<Type>::burstCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(p, iF),
    burstFvPatchField(p)
{}


template<class Type>
Foam::burstCyclicFvPatchField<Type>::burstCyclicFvPatchField
(
    const burstCyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    cyclicFvPatchField<Type>(ptf, p, iF, mapper),
    burstFvPatchField(ptf, p)
{}


template<class Type>
Foam::burstCyclicFvPatchField<Type>::burstCyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    cyclicFvPatchField<Type>(p, iF, dict),
    burstFvPatchField(p)
{}


template<class Type>
Foam::burstCyclicFvPatchField<Type>::burstCyclicFvPatchField
(
    const burstCyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    cyclicFvPatchField<Type>(ptf, iF),
    burstFvPatchField(ptf, ptf.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::burstCyclicFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    cyclicFvPatchField<Type>::autoMap(m);
    burstFvPatchField::autoMap(m);

}


template<class Type>
void Foam::burstCyclicFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    cyclicFvPatchField<Type>::rmap(ptf, addr);
    burstFvPatchField::rmap
    (
        dynamicCast<const burstCyclicFvPatchField<Type>>(ptf),
        addr
    );
}


template<class Type>
void Foam::burstCyclicFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    cyclicFvPatchField<Type>::updateCoeffs();
    burstFvPatchField::update();

}


template<class Type>
void Foam::burstCyclicFvPatchField<Type>::write(Ostream& os) const
{
    cyclicFvPatchField<Type>::write(os);
}


// ************************************************************************* //
