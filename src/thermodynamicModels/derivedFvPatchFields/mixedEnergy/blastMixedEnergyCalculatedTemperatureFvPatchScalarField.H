/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Changed class name
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
    Foam::blastMixedEnergyCalculatedTemperatureFvPatchScalarField

Description
    Base class for temperature boundary conditions in which the parameters of
    the mixed energy condition can be set directly.

SourceFiles
    blastMixedEnergyCalculatedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef blastMixedEnergyCalculatedTemperatureFvPatchScalarField_H
#define blastMixedEnergyCalculatedTemperatureFvPatchScalarField_H

#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class blastMixedEnergyCalculatedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class blastMixedEnergyCalculatedTemperatureFvPatchScalarField
:
    public calculatedFvPatchScalarField
{
    // Private Data

        //- Energy value field
        scalarField heRefValue_;

        //- Energy normal gradient field
        scalarField heRefGrad_;

        //- Fraction (0-1) of energy value used for boundary condition
        scalarField heValueFraction_;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("blastMixedEnergyCalculatedTemperature");


    // Constructors

        //- Construct from patch and internal field
        blastMixedEnergyCalculatedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blastMixedEnergyCalculatedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        blastMixedEnergyCalculatedTemperatureFvPatchScalarField
        (
            const blastMixedEnergyCalculatedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        blastMixedEnergyCalculatedTemperatureFvPatchScalarField
        (
            const blastMixedEnergyCalculatedTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new blastMixedEnergyCalculatedTemperatureFvPatchScalarField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        blastMixedEnergyCalculatedTemperatureFvPatchScalarField
        (
            const blastMixedEnergyCalculatedTemperatureFvPatchScalarField&,
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
                new blastMixedEnergyCalculatedTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Access

            //- Return true if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Access the energy value
            inline scalarField& heRefValue()
            {
                return heRefValue_;
            }

            //- Access the energy value
            inline const scalarField& heRefValue() const
            {
                return heRefValue_;
            }

            //- Access the energy gradient
            inline scalarField& heRefGrad()
            {
                return heRefGrad_;
            }

            //- Access the energy gradient
            inline const scalarField& heRefGrad() const
            {
                return heRefGrad_;
            }

            //- Access the energy fraction
            inline scalarField& heValueFraction()
            {
                return heValueFraction_;
            }

            //- Access the energy fraction
            inline const scalarField& heValueFraction() const
            {
                return heValueFraction_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            //  Used to update fields following mesh topology change
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            //  Used to reconstruct fields
            virtual void rmap(const fvPatchScalarField&, const labelList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
