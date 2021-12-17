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
    Foam::wedgeTractionFvPatchVectorField

Description
    Specifies contact traction for a wedge boundary patch according
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
        type            wedgeTraction;
        traction        (10 10 10);
        rampEndTime     1.0;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    wedgeTractionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeTractionFvPatchVectorField_H
#define wedgeTractionFvPatchVectorField_H

#include "wedgeFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class wedgeTractionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class wedgeTractionFvPatchVectorField
:
    public wedgeFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("wedgeTraction");


    // Constructors

        //- Construct from patch and internal field
        wedgeTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        wedgeTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given wedgeTractionFvPatchVectorField
        //  onto a new patch
        wedgeTractionFvPatchVectorField
        (
            const wedgeTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        wedgeTractionFvPatchVectorField
        (
            const wedgeTractionFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        wedgeTractionFvPatchVectorField
        (
            const wedgeTractionFvPatchVectorField&,
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
                new wedgeTractionFvPatchVectorField(*this, iF)
            );
        }

    // Member functions

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        // Write
        virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
