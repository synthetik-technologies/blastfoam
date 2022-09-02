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
    fixedDisplacementZeroShearFvPatchVectorField

Description
    The component of the displacement normal to the patch is applied,
    and the shear traction is set to zero.

    Patch contains non-orthogonal correction.

    Optionally time-varying by specifying a time series.

SourceFiles
    fixedDisplacementZeroShearFvPatchVectorField.C

Author
    Philip Cardiff, UCD/UT. All rights reserved.
    Based on velocityDisplacement by Aleksandar Karac

\*---------------------------------------------------------------------------*/

#ifndef fixedDisplacementZeroShearFvPatchVectorField_H
#define fixedDisplacementZeroShearFvPatchVectorField_H

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

class fixedDisplacementZeroShearFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Constant total displacement field
        vectorField totalDisp_;

        //- Displacement time series
        autoPtr<Function1<vector>> dispSeries_;

        //- Switch to force the shear gradient to be zero
        const Switch forceZeroShearGrad_;

public:

    //- Runtime type information
    TypeName("fixedDisplacementZeroShear");


    // Constructors

        //- Construct from patch and internal field
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        // fixedDisplacementZeroShearFvPatchVectorField
        //  onto a new patch
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        fixedDisplacementZeroShearFvPatchVectorField
        (
            const fixedDisplacementZeroShearFvPatchVectorField&,
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
                new fixedDisplacementZeroShearFvPatchVectorField(*this, iF)
            );
        }

    // Destructor

        virtual ~fixedDisplacementZeroShearFvPatchVectorField()
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
