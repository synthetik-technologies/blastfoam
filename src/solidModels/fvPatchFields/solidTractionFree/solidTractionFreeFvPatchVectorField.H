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

Author
   Philip Cardiff, UCD. All right reserved.

SourceFiles
    solidTractionFreeFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionFreeFvPatchVectorField_H
#define solidTractionFreeFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidTractionFreeFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionFreeFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
public:

    //- Runtime type information
    TypeName("solidTractionFree");


    // Constructors

        //- Construct from patch and internal field
        solidTractionFreeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionFreeFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidTractionFreeFvPatchVectorField onto a new patch
        solidTractionFreeFvPatchVectorField
        (
            const solidTractionFreeFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidTractionFreeFvPatchVectorField
        (
            const solidTractionFreeFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidTractionFreeFvPatchVectorField
        (
            const solidTractionFreeFvPatchVectorField&,
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
                new solidTractionFreeFvPatchVectorField(*this, iF)
            );
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
