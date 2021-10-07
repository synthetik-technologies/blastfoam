/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    displacementOrTractionFvPatchVectorField

Description
    This boundary condition allows a displacement or traction to specified on
    each face, where the valueFraction controls which is used.
    In addition, the valueFraction allows individual components to be specified.
    For example, it is possible to specify a normal displacement with a zero
    shear traction (like fixedDisplacementZeroShear), or a normal traction with
    a non-zero or zero shear displacement, or some faces could be traction and
    some could be displacement.

    It essentially a general solid boundary condition that can reproduce the
    behaviour of all other standard solid boundary conditions.

    The specifyNormalDirection field allows the normal direction to be used and
    takes precedence over the valueFraction field, using the following rules:
        1 means fix normal displacement
        0 means ignore this flag and use the valueFraction instead
        -1 means fix normal traction

    Patch contains non-orthogonal correction.

    Optionally time-varying by specifying a time series.

SourceFiles
    displacementOrTractionFvPatchVectorField.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef displacementOrTractionFvPatchVectorField_H
#define displacementOrTractionFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "Switch.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class directionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class displacementOrTractionFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Constant displacement field
        vectorField constantDisplacement_;

        //- Constant traction field
        vectorField constantTraction_;

        //- Optional: displacement time series
        autoPtr<Function1<vector>> displacementSeries_;

        //- Optional: traction time series
        autoPtr<Function1<vector>> tractionSeries_;

        //- Optional: fix displacement in normal or tangential directions
        //  1 means fix normal displacement
        //  0 means ignore this flag and use the valueFraction instead
        //  -1 means fix normal traction
        Field<label> specifyNormalDirection_;

public:

    //- Runtime type information
    TypeName("displacementOrTraction");


    // Constructors

        //- Construct from patch and internal field
        displacementOrTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        displacementOrTractionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // displacementOrTractionFvPatchVectorField
        //  onto a new patch
        displacementOrTractionFvPatchVectorField
        (
            const displacementOrTractionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        displacementOrTractionFvPatchVectorField
        (
            const displacementOrTractionFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        displacementOrTractionFvPatchVectorField
        (
            const displacementOrTractionFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<vector> > clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<vector> >
            (
                new displacementOrTractionFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~displacementOrTractionFvPatchVectorField()
        {}


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<vector>&,
                const labelList&
            );

        // Access functions

            // Const access to displacement series
            const Function1<vector>& displacementSeries() const
            {
                return displacementSeries_();
            }

            // Const access to tractionlacement series
            const Function1<vector>& tractionSeries() const
            {
                return tractionSeries_();
            }

            // Const access to specifyNormalDirection field
            const Field<label>& specifyNormalDirection() const
            {
                return specifyNormalDirection_;
            }

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
