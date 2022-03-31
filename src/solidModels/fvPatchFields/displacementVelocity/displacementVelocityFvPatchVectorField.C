/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "displacementVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcMeshPhi.H"
#include "pointMesh.H"
#include "pointFields.H"
#include "fixedValuePointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

displacementVelocityFvPatchVectorField::displacementVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


displacementVelocityFvPatchVectorField::displacementVelocityFvPatchVectorField
(
    const displacementVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


displacementVelocityFvPatchVectorField::displacementVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


displacementVelocityFvPatchVectorField::displacementVelocityFvPatchVectorField
(
    const displacementVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void displacementVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


//     if (internalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement for
        // this time-step
        Field<vector>::operator=
        (
            db().lookupObject<volVectorField>
            (
                "DD"
            ).boundaryField()[patch().index()]
           /db().time().deltaTValue()
        );
    }
//     else
//     {
//         // Lookup the old time total displacement
//         const volVectorField& D = db().lookupObject<volVectorField>("D");
//         const volVectorField& Dold = D.oldTime();
//
//         // The new total displacement is equal to Dold plus the increment of
//         // displacement based on the current velocity and time-step
//         Field<vector>::operator=
//         (
//             (
//                 D.boundaryField()[patch().index()]
//               - Dold.boundaryField()[patch().index()]
//             )/db().time().deltaTValue()
//         );
//     }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void displacementVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    displacementVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
