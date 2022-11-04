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

#include "hWaveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "uniformDimensionedFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::hWaveTransmissiveFvPatchField<Type>::hWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::hWaveTransmissiveFvPatchField<Type>::hWaveTransmissiveFvPatchField
(
    const hWaveTransmissiveFvPatchField<Type>& hwtpsf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(hwtpsf, p, iF, mapper)
{}


template<class Type>
Foam::hWaveTransmissiveFvPatchField<Type>::hWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::hWaveTransmissiveFvPatchField<Type>::hWaveTransmissiveFvPatchField
(
    const hWaveTransmissiveFvPatchField<Type>& hwtpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(hwtpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::hWaveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    scalarField pphi(this->size(), 0.0);
    if
    (
        this->patch().boundaryMesh().mesh().template foundObject
        <
            surfaceScalarField
        >(this->phiName_)
    )
    {
        pphi =
            this->patch().template lookupPatchField<surfaceScalarField, scalar>
            (
                this->phiName_
            );
    }
    else
    {
        WarningInFunction
            << this->phiName_ << " was not found" << endl;
    }

    const uniformDimensionedVectorField& g =
        this->db().template lookupObject<uniformDimensionedVectorField>("g");

    const fvPatchScalarField& ph =
        this->patch().template lookupPatchField<volScalarField, scalar>("h");

    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the gravitational wavespeed.
    return
        pphi/this->patch().magSf()
      + sqrt(mag(g.value())*(ph.patchInternalField()));
}


template<class Type>
void Foam::hWaveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    advectiveFvPatchField<Type>::write(os);
}


// ************************************************************************* //
