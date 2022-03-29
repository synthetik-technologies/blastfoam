/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::holeVelocityCyclicFvPatchVectorField

Description

Usage

See also
    Foam::cyclicFvPatchField

SourceFiles
    burstCyclicAMIFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef holeVelocityCyclicFvPatchVectorField_H
#define holeVelocityCyclicFvPatchVectorField_H

#include "cyclicFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class holeVelocityCyclicFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class holeVelocityCyclicFvPatchVectorField
:
    public cyclicFvPatchVectorField
{
    // Private Data

        //- Size of a hole
        scalar L_;

        //- Name of the length field
        word LName_;

        //- Mask (0 is wall, 1 is hole)
        scalarField mask_;

        //- Wall patchField
        tmp<fvPatchVectorField> wallPatchField_;


public:

    //- Runtime type information
    TypeName("holeVelocityCyclic");


    // Constructors

        //- Construct from patch and internal field
        holeVelocityCyclicFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        holeVelocityCyclicFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given holeVelocityCyclicFvPatchVectorField
        //  onto a new patch
        holeVelocityCyclicFvPatchVectorField
        (
            const holeVelocityCyclicFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        holeVelocityCyclicFvPatchVectorField
        (
            const holeVelocityCyclicFvPatchVectorField&
        ) = delete;

        //- Copy constructor setting internal field reference
        holeVelocityCyclicFvPatchVectorField
        (
            const holeVelocityCyclicFvPatchVectorField&,
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
                new holeVelocityCyclicFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Return neighbour field of internal field
            virtual tmp<vectorField> patchNeighbourField() const;


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap(const fvPatchVectorField&, const labelList&);


        // Evaluation functions

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType
            );

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
