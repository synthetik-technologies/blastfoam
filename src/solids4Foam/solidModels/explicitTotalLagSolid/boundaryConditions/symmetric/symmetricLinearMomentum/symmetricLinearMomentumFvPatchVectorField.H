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
    Foam::symmetricLinearMomentumFvPatchVectorField

Description
    Specifies contact linear momentum for a symmetric boundary patch
    according to the Riemann solver. External boundary traction can
    be specified together with a linear ramp function.

Usage
    \table
        Property     | Description              | Required    | Default value
        traction     | boundary traction        | no          | vector(0,0,0)
        rampEndTime  | ramp time for traction   | no          | VSMALL
    \endtable

    Example of the boundary condition specification:
    \verbatim
    inlet
    {
        type            symmetricLinearMomentum;
    }
    \endverbatim

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    symmetricLinearMomentumFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetricLinearMomentumFvPatchVectorField_H
#define symmetricLinearMomentumFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "fixedValueFvsPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class symmetricLinearMomentumFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class symmetricLinearMomentumFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
public:

    // Runtime type information
    TypeName("symmetricLinearMomentum");


    // Constructors

        //- Construct from patch and internal field
        symmetricLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetricLinearMomentumFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  symmetricLinearMomentumFvPatchVectorField onto a new patch
        symmetricLinearMomentumFvPatchVectorField
        (
            const symmetricLinearMomentumFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        symmetricLinearMomentumFvPatchVectorField
        (
            const symmetricLinearMomentumFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        symmetricLinearMomentumFvPatchVectorField
        (
            const symmetricLinearMomentumFvPatchVectorField&,
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
                new symmetricLinearMomentumFvPatchVectorField(*this, iF)
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
