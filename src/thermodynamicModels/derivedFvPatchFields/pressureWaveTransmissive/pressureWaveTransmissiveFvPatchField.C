/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

#include "pressureWaveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "fluidThermoModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::pressureWaveTransmissiveFvPatchField<Type>::
pressureWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    phaseName_(word::null)
{}


template<class Type>
Foam::pressureWaveTransmissiveFvPatchField<Type>::
pressureWaveTransmissiveFvPatchField
(
    const pressureWaveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    phaseName_(ptf.phaseName_)
{}


template<class Type>
Foam::pressureWaveTransmissiveFvPatchField<Type>::
pressureWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    phaseName_(dict.lookupOrDefault<word>("phaseName", word::null))
{}


template<class Type>
Foam::pressureWaveTransmissiveFvPatchField<Type>::
pressureWaveTransmissiveFvPatchField
(
    const pressureWaveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    phaseName_(ptpsf.phaseName_)
{}


template<class Type>
Foam::pressureWaveTransmissiveFvPatchField<Type>::
pressureWaveTransmissiveFvPatchField
(
    const pressureWaveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    phaseName_(ptpsf.phaseName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::pressureWaveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    const surfaceScalarField& phi =
        this->db().template lookupObject<surfaceScalarField>(this->phiName_);

    fvsPatchField<scalar> phip =
        this->patch().template
            lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop =
            this->patch().template
                lookupPatchField<volScalarField, scalar>(this->rhoName_);

        phip /= rhop;
    }

   const fluidThermoModel& thermo =
        this->db().template lookupObject<fluidThermoModel>
        (
            IOobject::groupName("basicThermo", phaseName_)
        );
    // Lookup the velocity and compressibility of the patch
    tmp<scalarField> cp
    (
        thermo.speedOfSound(this->patch().index())
    );
    // Calculate the speed of the field wave w
    // by summing the component of the velocity normal to the boundary
    // and the speed of sound (sqrt(cp)).
    return phip/this->patch().magSf() + cp;
}


template<class Type>
void Foam::pressureWaveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    writeEntryIfDifferent<word>(os, "phi", "phi", this->phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", this->rhoName_);
    writeEntryIfDifferent<word>(os, "phaseName", word::null, phaseName_);

    if (this->lInf_ > small)
    {
        writeEntry(os, "fieldInf", this->fieldInf_);
        writeEntry(os, "lInf", this->lInf_);
    }

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
