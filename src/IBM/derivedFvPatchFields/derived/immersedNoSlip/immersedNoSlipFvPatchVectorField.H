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

Class
    Foam::immersedNoSlipFvPatchVectorField

Description


SourceFiles
    immersedNoSlipFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef immersedNoSlipFvPatchVectorField_H
#define immersedNoSlipFvPatchVectorField_H

#include "immersedFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class immersedNoSlipFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class immersedNoSlipFvPatchVectorField
:
    public immersedFixedValueFvPatchField<vector>
{
public:

    //- Runtime type information
    TypeName("immersedNoSlip");

    // Constructors

        //- Construct from patch and internal field
        immersedNoSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        immersedNoSlipFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given
        // immersedNoSlipFvPatchVectorField onto a new patch
        immersedNoSlipFvPatchVectorField
        (
            const immersedNoSlipFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        immersedNoSlipFvPatchVectorField
        (
            const immersedNoSlipFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        immersedNoSlipFvPatchVectorField
        (
            const immersedNoSlipFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector>> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector>>
            (
                new immersedNoSlipFvPatchVectorField(*this, iF)
            );
        }

    // Member Functions

        // Evaluation functions

            //- Update the boundary values and set the patchInternal field
            virtual void updateCoeffs();

            //- Set the internal values
            virtual void setInternal();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
