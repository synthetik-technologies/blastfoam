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

#include "zeroGradientImmersedBoundaryFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroGradientImmersedBoundaryFvPatchField<Type>::zeroGradientImmersedBoundaryFvPatchField
(
    GeometricField<Type, fvPatchField, volMesh>& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
:
    immersedBoundaryFvPatchField<Type>(f, dict, ibo),
    nSmooth_(dict.lookupOrDefault("nSmooth", 0))
{
    this->values_.resize(0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void
Foam::zeroGradientImmersedBoundaryFvPatchField<Type>::addForcing
(
    Field<Type>& F,
    const Field<scalar>& alphaRho,
    const Field<Type>& old,
    const Field<Type>& RHS,
    const scalar& dt
) const
{}


template<class Type>
void Foam::zeroGradientImmersedBoundaryFvPatchField<Type>::updateCoeffs() const
{
    const labelList& pE(this->ibm_.patchExternalCells());
    forAll(pE, i)
    {
        const label celli = pE[i];
        if (celli >= 0)
        {
            this->values_[i] = this->field_[celli];
        }
    }
}


template<class Type>
void Foam::zeroGradientImmersedBoundaryFvPatchField<Type>::setValues()
{
    const labelList& pI(this->ibm_.patchInternalCells());
    const labelList& pE(this->ibm_.patchExternalCells());
    const scalarList& magSf(this->ibm_.magSf());
    Field<Type> newValues(this->ibm_.pMesh().nCells(), Zero);

    scalarField w(newValues.size(), 0);
    forAll(pI, i)
    {
        const label celli = pI[i];
        const label cellj = pE[i];
        if (celli >= 0 && cellj >= 0)
        {
            w[celli] += magSf[i];
            newValues[celli] = this->field_[cellj];//*magSf[i];
        }
    }
    for (label i = 0; i < nSmooth_; i++)
    {
        this->setSmoothInternal();
    }

    forAll(pI, i)
    {
        const label celli = pI[i];
        const label cellj = pE[i];
        if (celli >= 0 && cellj >= 0)
        {
            this->field_[celli] = newValues[celli];///w[celli];
        }
    }
}

// ************************************************************************* //
