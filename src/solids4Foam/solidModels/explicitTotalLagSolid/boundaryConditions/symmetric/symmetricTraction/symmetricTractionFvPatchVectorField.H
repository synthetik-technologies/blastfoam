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
    Foam::symmetricTractionFvPatchVectorField

Description
    Specifies contact traction for a symmetric boundary patch according
    to the Riemann solver. External boundary traction can also be specified.

Usage
    \table
        Property    | Description               | Required  | Default value
        traction    | boundary traction         | no        | vector(0,0,0)
        rampEndTime | ramp time for traction    | no        | VSMALL

    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            symmetricTraction;
        traction        (10 10 10);
        rampEndTime     1.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    symmetricTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetricTractionFvPatchVectorField_H
#define symmetricTractionFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class symmetricTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class symmetricTractionFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("symmetricTraction");


    // Constructors

        //- Construct from patch and internal field
        symmetricTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetricTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given symmetricTractionFvPatchVectorField
        //  onto a new patch
        symmetricTractionFvPatchVectorField
        (
            const symmetricTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        symmetricTractionFvPatchVectorField
        (
            const symmetricTractionFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        symmetricTractionFvPatchVectorField
        (
            const symmetricTractionFvPatchVectorField&,
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
                new symmetricTractionFvPatchVectorField(*this, iF)
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
