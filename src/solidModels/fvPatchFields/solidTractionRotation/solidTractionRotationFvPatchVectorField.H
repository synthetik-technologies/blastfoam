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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    solidTractionRotationFvPatchVectorField

Description
    Traction free boundary condition for solid foam solvers.

    solidTraction boundary condition where the traction and pressure are set to
    zero.

Author
   Philip Cardiff, UCD. All right reserved.

SourceFiles
    solidTractionRotationFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef solidTractionRotationFvPatchVectorField_H
#define solidTractionRotationFvPatchVectorField_H

#include "fvPatchFields.H"
#include "solidTractionFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
             Class solidTractionRotationFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class solidTractionRotationFvPatchVectorField
:
    public solidTractionFvPatchVectorField
{
    //- Origin or rotation
    vector rotationOrigin_;

    //- Axis of rotation
    vector rotationAxis_;

    //- Original face Centres
    pointField origFaceCentres_;

    //- Original patch points
    pointField origPatchPoints_;


public:

    //- Runtime type information
    TypeName("solidTractionRotation");


    // Constructors

        //- Construct from patch and internal field
        solidTractionRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        solidTractionRotationFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  solidTractionRotationFvPatchVectorField onto a new patch
        solidTractionRotationFvPatchVectorField
        (
            const solidTractionRotationFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        solidTractionRotationFvPatchVectorField
        (
            const solidTractionRotationFvPatchVectorField&
        ) = delete;

        //- Construct as copy setting internal field reference
        solidTractionRotationFvPatchVectorField
        (
            const solidTractionRotationFvPatchVectorField&,
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
                new solidTractionRotationFvPatchVectorField(*this, iF)
            );
        }

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );

            virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
