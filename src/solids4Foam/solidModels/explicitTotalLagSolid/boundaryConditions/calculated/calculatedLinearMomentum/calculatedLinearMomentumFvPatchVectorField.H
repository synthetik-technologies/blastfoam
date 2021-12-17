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
    Foam::calculatedLinearMomentumFvPatchVectorField

Description
    Linear momentum is calculated from U and rho

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    calculatedLinearMomentumPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedLinearMomentumFvPatchVectorField_H
#define calculatedLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class calculatedLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class calculatedLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("calculatedLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        calculatedLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        calculatedLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  calculatedLinearMomentumFvPatchVectorField onto a new patch
        calculatedLinearMomentumFvPatchVectorField
        (
            const calculatedLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        calculatedLinearMomentumFvPatchVectorField
        (
            const calculatedLinearMomentumFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        calculatedLinearMomentumFvPatchVectorField
        (
            const calculatedLinearMomentumFvPatchVectorField&,
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
                new calculatedLinearMomentumFvPatchVectorField
                (
                    *this,
                    iF
                )
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
