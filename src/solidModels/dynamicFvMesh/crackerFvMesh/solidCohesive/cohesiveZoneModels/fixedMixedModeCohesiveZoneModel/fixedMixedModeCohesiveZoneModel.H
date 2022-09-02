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
    fixedMixedModeCohesiveZoneModel

Description
    Mixed-Mode cohesiv zone model as detailed in:
    Camacho GT, Ortiz M (1996) Computational modelling of impact damage in
    brittle materials. Int J Solids Struct 33:2899–2938.

    This model assumes the mode-mixity of a damage face to remain constant
    from initiation to final failure.

    Linear softening cohesive zone shape is assumed; Dugdale can be specified
    as an option.

Author
    Zeljko Tukovic, FSB Zagreb. All rights reserved.
    Re-organised by Philip Cardif, UCD/UT. All rights reserved.

SourceFiles
    fixedMixedModeCohesiveZoneModel.C

\*---------------------------------------------------------------------------*/

#ifndef fixedMixedModeCohesiveZoneModel_H
#define fixedMixedModeCohesiveZoneModel_H

#include "cohesiveZoneModel.H"
#include "solidCohesiveFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fixedMixedModeCohesiveZoneModel Declaration
\*---------------------------------------------------------------------------*/

class fixedMixedModeCohesiveZoneModel
:
    public cohesiveZoneModel
{
    // Private data

        //- Strength in mode-I crack/damage
        const dimensionedScalar sigmaMax_;

        //- Critical mode-I fracture energy crack/damage
        const dimensionedScalar GIc_;

        //- Critical opening delta
        //  Note:
        //  linear shape: deltaC = 2.0*GIc/sigmaMax
        //  Dugdale shape: deltaC = GIc/sigmaMax
        dimensionedScalar deltaC_;

        //- If a face has cracked:
        //  1 = cracked
        //  0 = damaged but not cracked
        Field<bool> cracked_;

        //- Initiation traction: this is the traction the face initiated at
        vectorField initTraction_;

        //- Current delta normal component
        scalarField deltaN_;

        //- Current delta shear component
        vectorField deltaS_;

        //- Current delta effective magnitude
        scalarField deltaEff_;

        //- Current unloading separation distance
        vectorField unloadingDeltaEff_;

        //- Minimal unloading relative separation distance: this is a fraction
        //  of the critical delta
        scalar minUnloadingDelta_;

        //- Ratio between shear and normal critical traction
        const scalar beta_;

        //- Switch for cohesive zone model shape: dugdale or linear softening
        const Switch dugdale_;

        //- Penalty factor scale
        const scalar penaltyScale_;

        //- Penalty factor for penalty contact
        mutable scalar penaltyFactor_;


    // Private Member Functions

        // Calculate traction due to contact for the current face
        vector contactTraction
        (
            const vector& delta,
            const scalar deltaN,
            const vector& faceN,
            const scalar faceFricCoeff
        ) const;

        // Calculate traction in the damage zone during loading
        vector damageTraction
        (
            const scalar faceDeltaEff,
            const scalar faceSigmaMax,
            const scalar faceDeltaC,
            const vector& faceInitTrac
        ) const;

        // Calculate traction in the damage zone during unloading
        vector unloadingDamageTraction
        (
            const scalar faceDeltaEff,
            const scalar faceSigmaMax,
            const scalar faceDeltaC,
            const vector& faceInitTrac,
            const scalar faceMagUnloadingDeltaEff
        ) const;

        //- Calculate and set the penalty factor
        void calcPenaltyFactor() const;

        //- Disallow default bitwise assignment
        void operator=(const fixedMixedModeCohesiveZoneModel&);


public:

    //- Runtime type information
    TypeName("fixedMixedMode");

    // Static data members


    // Constructors

        //- Construct from dictionary
        fixedMixedModeCohesiveZoneModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict
        );

        //- Construct as copy
        fixedMixedModeCohesiveZoneModel(const fixedMixedModeCohesiveZoneModel&);

        //- Construct and return a clone
        virtual autoPtr<cohesiveZoneModel> clone() const
        {
            return
                autoPtr<cohesiveZoneModel>
                (
                    new fixedMixedModeCohesiveZoneModel(*this)
                );
        }


    // Destructor

        virtual ~fixedMixedModeCohesiveZoneModel();


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper& mapper
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const solidCohesiveFvPatchVectorField& sc,
                const labelList& addr
            );

        // Access functions

            //- Update old values
            virtual void updateOldFields();

            //- Update fracture energy
            virtual void updateEnergy
            (
                const vectorField& traction,
                const vectorField& delta
            );

            //- Update cohesive traction field
            virtual void updateTraction
            (
                vectorField& traction,
                const vectorField& delta
            );

            //- Return traction fraction for undamaged faces
            //  This function is called by the faceBreakerLaw to determine which
            //  internal faces should be add to the cohesive zone
            virtual tmp<surfaceScalarField> initiationTractionFraction() const;

            //- crack and damage field for visualisation
            //- cracked faces set to 2, damaged but not cracked set to 1
            virtual tmp<scalarField> crackingAndDamage() const;

            //- Write dictionary
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
