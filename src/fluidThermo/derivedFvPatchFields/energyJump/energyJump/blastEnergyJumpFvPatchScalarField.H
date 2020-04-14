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

Class
    Foam::blastEnergyJumpFvPatchScalarField

Description
    This boundary condition provides an energy jump condition across a pair
    of coupled patches.  It is not applied directly, but is employed on-the-fly
    when converting temperature boundary conditions into energy.

See also
    Foam::fixedJumpFvPatchField

SourceFiles
    blastEnergyJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef blastEnergyJumpFvPatchScalarField_H
#define blastEnergyJumpFvPatchScalarField_H

#include "fixedJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class blastEnergyJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class blastEnergyJumpFvPatchScalarField
:
    public fixedJumpFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("blastEnergyJump");

    // Constructors

        //- Construct from patch and internal field
        blastEnergyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blastEnergyJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blastEnergyJumpFvPatchScalarField onto a
        //  new patch
        blastEnergyJumpFvPatchScalarField
        (
            const blastEnergyJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        blastEnergyJumpFvPatchScalarField
        (
            const blastEnergyJumpFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new blastEnergyJumpFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        blastEnergyJumpFvPatchScalarField
        (
            const blastEnergyJumpFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar>> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar>>
            (
                new blastEnergyJumpFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
