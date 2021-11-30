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
    simpleCohesiveZoneFvPatchVectorField

Description
    Cohesive zone boundary condition to be used in conjunction with the
    simpleCrackerFvMesh dynamicFvMesh class.

    This boundary condition is for mode-I fracture only, where the fracture
    propagates along a symmetryPlane.

    Essentially, this boundary condition is a symmetryPlane that transfer
    faces to a cohesiveZone when a certain traction criteria is met.

SourceFiles
    simpleCohesiveZoneFvPatchVectorField.C

Author
    Zeljko Tukovic. FSB Zaregb. All rights reserved.
    Philip Cardiff. UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef simpleCohesiveZoneFvPatchVectorField_H
#define simpleCohesiveZoneFvPatchVectorField_H

#ifdef FOAMEXTEND
    #include "foamTime.H"
#endif
#include "fvPatchFields.H"
#include "solidDirectionMixedFvPatchVectorField.H"
#include "simpleCohesiveZoneLaw.H"
#include "Switch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class solidDirectionMixedFvPatch Declaration
\*---------------------------------------------------------------------------*/

class simpleCohesiveZoneFvPatchVectorField
:
    public solidDirectionMixedFvPatchVectorField
{
    // Private data

        //- Copy of the input dictionary
        const dictionary dict_;

        //- Accumulated displacement
        vectorField totRefValue_;

        //- Cohesive law
        mutable autoPtr<simpleCohesiveZoneLaw> cohesiveLawPtr_;

        //- Crack indicator
        scalarField crackIndicator_;

        //- Damage/craze indicator
        scalarField damageIndicator_;

        //- Under-relaxation factor
        scalar relaxationFactor_;

        //- Separation distance
        scalarField separationDistance_;

        //- Olde separation distance
        scalarField oldSeparationDistance_;

        //- Current unloading separation distance
        scalarField unloadingSeparationDistance_;

        //- Use explicit spearation distance
        Switch explicitSeparationDistance_;

        //- Current time index
        label curTimeIndex_;

        //- Initiation traction
        vectorField initiationTraction_;

        //- Break only one face per updateCrack
        Switch breakOnlyOneFace_;

    //- Private Member Functions

        //- Make the cohesive law
        void makeSimpleCohesiveZoneLaw() const;

        //- Return the cohesive law
        const simpleCohesiveZoneLaw& law() const;

    public:

    //- Runtime type information
    TypeName("simpleCohesiveZone");


    // Constructors

        //- Construct from patch and internal field
        simpleCohesiveZoneFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        simpleCohesiveZoneFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  simpleCohesiveZoneFvPatchVectorField onto
        //  a new patch
        simpleCohesiveZoneFvPatchVectorField
        (
            const simpleCohesiveZoneFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy setting internal field reference
        simpleCohesiveZoneFvPatchVectorField
        (
            const simpleCohesiveZoneFvPatchVectorField&,
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
                new simpleCohesiveZoneFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return total displacement
            const vectorField& totRefValue() const
            {
                return totRefValue_;
            }

            //- Return crack indicator
            const scalarField& crackIndicator() const
            {
                return crackIndicator_;
            }

            //- Return DAMAGE indicator
            const scalarField& damageIndicator() const
            {
                return damageIndicator_;
            }

            //- Return separation distance
            const scalarField& separationDistance() const
            {
                return separationDistance_;
            }

            //- Return relative separation distance
            const tmp<scalarField> relSeparationDistance() const
            {
                return tmp<scalarField>
                (
                    new scalarField(separationDistance_/law().deltaC().value())
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
                const fvPatchField<vector>&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Update crack
            label updateCrack();

//             //- Evaluate the patch field
//             virtual void evaluate
//             (
//                 const Pstream::commsTypes commsType=Pstream::blocking
//             );

        //- Write
        virtual void write(Ostream&) const;


    // Member operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
