/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Class
    Foam::calculatedTractionFvPatchVectorField

Description
    Specifies contact traction for a traction boundary patch according to
    the Riemann solver. External boundary traction can also be specified.

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    calculatedTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedTractionFvPatchVectorField_H
#define calculatedTractionFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fvsPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class calculatedTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class calculatedTractionFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("calculatedTraction");


    // Constructors

        //- Construct from patch and internal field
        calculatedTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        calculatedTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given calculatedTractionFvPatchVectorField
        //  onto a new patch
        calculatedTractionFvPatchVectorField
        (
            const calculatedTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        calculatedTractionFvPatchVectorField
        (
            const calculatedTractionFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        calculatedTractionFvPatchVectorField
        (
            const calculatedTractionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new calculatedTractionFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        // Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        // Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
