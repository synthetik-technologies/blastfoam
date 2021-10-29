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

#include "slipImmersedBoundaryFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::slipImmersedBoundaryFvPatchVectorField::
slipImmersedBoundaryFvPatchVectorField
(
    volVectorField& f,
    const dictionary& dict,
    const immersedBoundaryObject& ibo
)
:
    immersedBoundaryVectorPatchField(f, dict, ibo)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::slipImmersedBoundaryFvPatchVectorField::addForcing
(
    Field<vector>& F,
    const Field<scalar>& alphaRho,
    const Field<vector>& old,
    const Field<vector>& RHS,
    const scalar& dt
) const
{
    //- Only add forcing in the normal direction
    vectorField normal(this->ibm_.Sf()/this->ibm_.magSf());
    vectorField wallU(ibm_.velocity(ibm_.faceCentres()));
    tmp<Field<vector>> interpF
    (
        (
            wallU*ibm_.interpolateTo(alphaRho)
          - ibm_.interpolateTo(old)
        )/dt
      + ibm_.interpolateTo(RHS)

    );
    interpF.ref() = (interpF() & normal)*normal;
    ibm_.interpolateFrom(interpF(), F);
}


void Foam::slipImmersedBoundaryFvPatchVectorField::updateCoeffs() const
{
    vectorField v(this->ibm_.interpolateTo(this->field_));
    vectorField normal(this->ibm_.Sf()/this->ibm_.magSf());
    vectorField wallU(ibm_.velocity(ibm_.faceCentres()));

    tmp<vectorField> vt(v - (v & normal)*normal);
    tmp<vectorField> vn((wallU & normal)*normal);

    values_ = vn + vt;
}


void Foam::slipImmersedBoundaryFvPatchVectorField::setValues()
{
//     this->ibm_.setInternal
//     (
//         field_,
//         this->ibm_.velocity(ibm_.internalC())()
//     );

    const labelList& pI(this->ibm_.patchInternalCells());
    const labelList& pE(this->ibm_.patchExternalCells());
    labelList n(pE.size(), 0);
    forAll(pI, i)
    {
        const label celli = pI[i];
        if (celli >= 0)
        {
            if (n[i]++ == 0)
            {
                this->field_[celli] = this->values_[i];
            }
            else
            {
                this->field_[celli] += this->values_[i];
            }
        }
    }
    forAll(pE, i)
    {
        if (n[i] > 0)
        {
            this->field_[pI[i]] /= scalar(n[i]);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeImmersedPatchTypeField
    (
        immersedBoundaryVectorPatchField,
        slipImmersedBoundaryFvPatchVectorField
    );
}


// ************************************************************************* //
