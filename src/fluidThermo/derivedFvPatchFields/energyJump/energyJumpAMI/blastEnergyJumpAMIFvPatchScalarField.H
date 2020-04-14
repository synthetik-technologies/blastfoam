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

See also
    Foam::fixedJumpAMIFvPatchField

SourceFiles
    blastEnergyJumpAMIFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef blastEnergyJumpAMIFvPatchScalarField_H
#define blastEnergyJumpAMIFvPatchScalarField_H

#include "fixedJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class blastEnergyJumpAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class blastEnergyJumpAMIFvPatchScalarField
:
    public fixedJumpAMIFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("blastEnergyJumpAMI");

    // Constructors

        //- Construct from patch and internal field
        blastEnergyJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blastEnergyJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blastEnergyJumpAMIFvPatchScalarField onto a
        //  new patch
        blastEnergyJumpAMIFvPatchScalarField
        (
            const blastEnergyJumpAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        blastEnergyJumpAMIFvPatchScalarField
        (
            const blastEnergyJumpAMIFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar>> clone() const
        {
            return tmp<fvPatchField<scalar>>
            (
                new blastEnergyJumpAMIFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        blastEnergyJumpAMIFvPatchScalarField
        (
            const blastEnergyJumpAMIFvPatchScalarField&,
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
                new blastEnergyJumpAMIFvPatchScalarField(*this, iF)
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
