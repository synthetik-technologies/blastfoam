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
    Foam::tractionLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a traction boundary patch according
    to the Riemann solver. External boundary traction or pressure can also
    be specified.

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    tractionLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef tractionLinearMomentumFvPatchVectorField_H
#define tractionLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class tractionLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class tractionLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("tractionLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        tractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given tractionLinearMomentumFvPatchVectorField
        //  onto a new patch
        tractionLinearMomentumFvPatchVectorField
        (
            const tractionLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        tractionLinearMomentumFvPatchVectorField
        (
            const tractionLinearMomentumFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        tractionLinearMomentumFvPatchVectorField
        (
            const tractionLinearMomentumFvPatchVectorField&,
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
                new tractionLinearMomentumFvPatchVectorField(*this, iF)
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
