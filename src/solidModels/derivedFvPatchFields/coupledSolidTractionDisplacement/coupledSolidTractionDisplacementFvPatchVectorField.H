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
    Foam::coupledSolidTractionDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation in which the traction is caused by
    the coupledSolidTraction pressure of an external liquid reservoir.

SourceFiles
    coupledSolidTractionDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef coupledSolidTractionDisplacementFvPatchVectorField_H
#define coupledSolidTractionDisplacementFvPatchVectorField_H

#include "tractionDisplacementFvPatchVectorField.H"
#include "mappedPatchSelector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class coupledSolidTractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledSolidTractionDisplacementFvPatchVectorField
:
    public tractionDisplacementFvPatchVectorField
{
    // Private Data

        //- Selector for mapped patch type
        mappedPatchSelector mpp_;

        //- Name of neighbour pressure field
        const word pName_;

        //- Reference pressure
        scalar pRef_;


public:

    //- Runtime type information
    TypeName("coupledSolidTraction");


    // Constructors

        //- Construct from patch and internal field
        coupledSolidTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coupledSolidTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  coupledSolidTractionDisplacementFvPatchVectorField onto a new patch
        coupledSolidTractionDisplacementFvPatchVectorField
        (
            const coupledSolidTractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        coupledSolidTractionDisplacementFvPatchVectorField
        (
            const coupledSolidTractionDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new coupledSolidTractionDisplacementFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        coupledSolidTractionDisplacementFvPatchVectorField
        (
            const coupledSolidTractionDisplacementFvPatchVectorField&,
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
                new coupledSolidTractionDisplacementFvPatchVectorField(*this, iF)
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
