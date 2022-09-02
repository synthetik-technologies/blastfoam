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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidTractionFreeFvPatchVectorField

Description
    Traction free boundary condition for solid foam solvers.

    solidTraction boundary condition where the traction and pressure are set to
    zero.

SourceFiles
    hydrostaticPressureFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef hydrostaticPressureFvPatchVectorField_H
#define hydrostaticPressureFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class hydrostaticPressureFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class hydrostaticPressureFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    //- Reference pressure
    scalar pRef_;

    //- Reference depth
    scalar hRef_;

    //- Reference density
    scalar rho_;

public:

    //- Runtime type information
    TypeName("hydrostaticPressure");


    // Constructors

        //- Construct from patch and internal field
        hydrostaticPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        hydrostaticPressureFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  hydrostaticPressureFvPatchVectorField onto a new patch
        hydrostaticPressureFvPatchVectorField
        (
            const hydrostaticPressureFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        hydrostaticPressureFvPatchVectorField
        (
            const hydrostaticPressureFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        hydrostaticPressureFvPatchVectorField
        (
            const hydrostaticPressureFvPatchVectorField&,
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
                new hydrostaticPressureFvPatchVectorField(*this, iF)
            );
        }

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
