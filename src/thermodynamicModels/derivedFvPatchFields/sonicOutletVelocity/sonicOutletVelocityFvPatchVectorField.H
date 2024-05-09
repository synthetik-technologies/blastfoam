/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::sonicOutletVelocityFvPatchVectorField

Description


See also
    Foam::thermoBasePatchField
    Foam::fixedValueFvPatchField

SourceFiles
    sonicOutletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef sonicOutletVelocityFvPatchVectorField_H
#define sonicOutletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "thermoBasePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class sonicOutletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class sonicOutletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField,
    public thermoBasePatchField
{
protected:

    // Protected Data

        //- Reference pressure
        const scalar pRef_;


public:

    //- Runtime type information
    TypeName("sonicOutletVelocity");


    // Constructors

        //- Construct from patch and internal field
        sonicOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        sonicOutletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given sonicOutletVelocityFvPatchVectorField
        //  onto a new patch
        sonicOutletVelocityFvPatchVectorField
        (
            const sonicOutletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Disallow copy without setting internal field reference
        sonicOutletVelocityFvPatchVectorField
        (
            const sonicOutletVelocityFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        sonicOutletVelocityFvPatchVectorField
        (
            const sonicOutletVelocityFvPatchVectorField&,
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
                new sonicOutletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
