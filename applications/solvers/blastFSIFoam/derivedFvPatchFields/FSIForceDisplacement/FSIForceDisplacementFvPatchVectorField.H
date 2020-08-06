/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2019 OpenFOAM Foundation
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
    Foam::FSIForceDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation in which the traction is caused by
    the FSIForce pressure of an external liquid reservoir.

SourceFiles
    FSIForceDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef FSIForceDisplacementFvPatchVectorField_H
#define FSIForceDisplacementFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class FSIForceDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class FSIForceDisplacementFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{
    // Private Data

        //- Name of neighbour pressure field
        const word pName_;

        //- Reference pressure
        scalarField pRef_;


public:

    //- Runtime type information
    TypeName("FSIForceDisplacement");


    // Constructors

        //- Construct from patch and internal field
        FSIForceDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        FSIForceDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  FSIForceDisplacementFvPatchVectorField onto a new patch
        FSIForceDisplacementFvPatchVectorField
        (
            const FSIForceDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        FSIForceDisplacementFvPatchVectorField
        (
            const FSIForceDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new FSIForceDisplacementFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        FSIForceDisplacementFvPatchVectorField
        (
            const FSIForceDisplacementFvPatchVectorField&,
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
                new FSIForceDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Map (and resize as needed) from self given a mapping object
        //  Used to update fields following mesh topology change
        virtual void autoMap(const fvPatchFieldMapper&);

        //- Reverse map the given fvPatchField onto this fvPatchField
        //  Used to reconstruct fields
        virtual void rmap(const fvPatchVectorField&, const labelList&);

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
