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

#include "subcriticalFixedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::subcriticalFixedFvPatchField<Type>::subcriticalFixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::subcriticalFixedFvPatchField<Type>::subcriticalFixedFvPatchField
(
    const subcriticalFixedFvPatchField<Type>& hwtpsf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(hwtpsf, p, iF, mapper)
{}


template<class Type>
Foam::subcriticalFixedFvPatchField<Type>::subcriticalFixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::subcriticalFixedFvPatchField<Type>::subcriticalFixedFvPatchField
(
    const subcriticalFixedFvPatchField<Type>& hwtpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(hwtpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::subcriticalFixedFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalarField pphi(this->size(), Zero);
    if
    (
        this->patch().boundaryMesh().mesh().template foundObject
        <
            surfaceScalarField
        >("phi")
    )
    {
        pphi =
            this->patch().template lookupPatchField<surfaceScalarField, scalar>
            (
                "phi"
            );
    }
    else
    {
        WarningInFunction
            << "phi was not found" << endl;
    }

    const fvPatchScalarField& ph =
        this->patch().template lookupPatchField<volScalarField, scalar>("h");
    tmp<scalarField> h =
        isA<mixedFvPatchField<scalar>>(ph)
      ? dynamicCast<const mixedFvPatchField<scalar>>(ph).refValue()
      : static_cast<const scalarField&>(ph);
    scalarField ws
    (
        max
        (
            sqrt
            (
                mag
                (
                    this->db().template
                        lookupObject<uniformDimensionedVectorField>("g")
                ).value()
               *h
            ),
            small
        )
    );

    scalarField Fr(mag(pphi)/(this->patch().magSf()*ws));

    this->valueFraction() = pos(1.0 - Fr);

    mixedFvPatchField<Type>::updateCoeffs();
}

// ************************************************************************* //
