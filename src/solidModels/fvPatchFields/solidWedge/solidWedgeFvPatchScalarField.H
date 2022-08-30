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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::solidWedgeFvPatchScalarField

Description
    Standard wedge boundary with non-orthogonal correction enabled.

SourceFiles
    solidWedgeFvPatchScalarField.C

Author
    Philip Cardiff, UCD. All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef solidWedgeFvPatchScalarField_H
#define solidWedgeFvPatchScalarField_H

#include "wedgeFvPatchFields.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class wedgeFvPatch Declaration
\*---------------------------------------------------------------------------*/

class solidWedgeFvPatchScalarField
:
    public wedgeFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("solidWedge");


    // Constructors

        //- Construct from patch and internal field
        solidWedgeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidWedgeFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidWedgeFvPatchScalarField onto a new
        //- patch
        solidWedgeFvPatchScalarField
        (
            const solidWedgeFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidWedgeFvPatchScalarField
        (
            const solidWedgeFvPatchScalarField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidWedgeFvPatchScalarField
        (
            const solidWedgeFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new solidWedgeFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<Field<scalar> > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


#endif

// ************************************************************************* //
