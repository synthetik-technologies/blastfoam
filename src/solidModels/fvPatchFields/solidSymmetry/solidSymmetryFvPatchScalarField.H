/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::solidSymmetryFvPatchScalarField

Description
    Symmetry plane with non-orthogonal correction for scalar field.

SourceFiles
    solidSymmetryFvPatchScalarField.C

Author
    Philip Cardiff, UCD. All rights reserved based on symmetryDisplacement

\*---------------------------------------------------------------------------*/

#ifndef solidSymmetryFvPatchScalarField_H
#define solidSymmetryFvPatchScalarField_H

#include "symmetryFvPatchFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class solidSymmetryFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class solidSymmetryFvPatchScalarField
:
    public symmetryFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("solidSymmetry");


    // Constructors

        //- Construct from patch and internal field
        solidSymmetryFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidSymmetryFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given solidSymmetryFvPatchScalarField
        //  onto a new patch
        solidSymmetryFvPatchScalarField
        (
            const solidSymmetryFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidSymmetryFvPatchScalarField
        (
            const solidSymmetryFvPatchScalarField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidSymmetryFvPatchScalarField
        (
            const solidSymmetryFvPatchScalarField&,
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
                new solidSymmetryFvPatchScalarField(*this, iF)
            );
        }

        tmp<Field<scalar> > snGrad() const;

        void evaluate(const Pstream::commsTypes);

        void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
