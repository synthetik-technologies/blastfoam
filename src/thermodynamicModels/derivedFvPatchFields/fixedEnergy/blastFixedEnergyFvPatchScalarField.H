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
    Foam::fixedEnergyFvPatchScalarField

Description
    This boundary condition provides a fixed condition for internal energy

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            blastFixedEnergy;
        value           uniform 100;
    }
    \endverbatim

See also
    Foam::blastFixedEnergyFvPatchScalarField

SourceFiles
    blastFixedEnergyFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef blastFixedEnergyFvPatchScalarField_H
#define blastFixedEnergyFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class blastFixedEnergyFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class blastFixedEnergyFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("blastFixedEnergy");


    // Constructors

        //- Construct from patch and internal field
        blastFixedEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blastFixedEnergyFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given blastFixedEnergyFvPatchScalarField
        // onto a new patch
        blastFixedEnergyFvPatchScalarField
        (
            const blastFixedEnergyFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        blastFixedEnergyFvPatchScalarField
        (
            const blastFixedEnergyFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new blastFixedEnergyFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        blastFixedEnergyFvPatchScalarField
        (
            const blastFixedEnergyFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new blastFixedEnergyFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
