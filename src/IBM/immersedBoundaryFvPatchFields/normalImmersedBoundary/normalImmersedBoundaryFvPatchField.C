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

#include "immersedBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::immersedBoundaryFvPatchField<Type>::immersedBoundaryFvPatchField
(
    GeometricField<Type, fvPatchField, volMesh>& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
:
    ibm_(ibo),
    field_(f),
    values_(ibo.nFaces(), Zero)
{}


template<class Type>
Foam::autoPtr<Foam::immersedBoundaryFvPatchField<Type>> Foam::immersedBoundaryFvPatchField<Type>::New
(
    GeometricField<Type, fvPatchField, volMesh>& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
{
    const word patchType(dict.lookup("type"));
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(patchType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown patchField type "
            << patchType << nl << nl
            << "Valid patchField types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(f, dict, ibo);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void
Foam::immersedBoundaryFvPatchField<Type>::addForcing
(
    Field<Type>& F,
    const Field<scalar>& alphaRho,
    const Field<Type>& old,
    const Field<Type>& RHS,
    const scalar& dt
) const
{
    updateCoeffs();

    tmp<Field<Type>> interpF
    (
        (
            values_*ibm_.interpolateTo(alphaRho)
          - ibm_.interpolateTo(old)
        )/dt
      + ibm_.interpolateTo(RHS)

    );
    ibm_.interpolateFrom(interpF(), F);
}


// ************************************************************************* //
