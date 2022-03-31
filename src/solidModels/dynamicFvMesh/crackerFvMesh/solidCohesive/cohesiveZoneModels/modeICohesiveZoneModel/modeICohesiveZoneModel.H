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
    modeICohesiveZoneModel

Description
    Mode-I Dugdale-type cohesive zone model.

Author
    Philip Cardif, UCD/UT. All rights reserved.

SourceFiles
    modeICohesiveZoneModel.C

\*---------------------------------------------------------------------------*/

#ifndef modeICohesiveZoneModel_H
#define modeICohesiveZoneModel_H

#include "cohesiveZoneModel.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Forward declaration
    class solidCohesiveFvPatchVectorField;

/*---------------------------------------------------------------------------* \
                         Class modeICohesiveZoneModel Declaration
\*---------------------------------------------------------------------------*/

class modeICohesiveZoneModel
:
    public cohesiveZoneModel
{
    // Private data

        //- Strength in mode-I crack/damage
        dimensionedScalar sigmaMax_;

        //- Critical mode-I fracture energy crack/damage
        dimensionedScalar GIc_;

        //- If a face has cracked:
        //  PC: it might make sense to move this field to the base class
        //  1 = cracked
        //  0 = damaged but not cracked
        Field<bool> cracked_;

        //- Current patch normal traction
        scalarField tractionN_;

        //- Old time patch normal traction
        scalarField oldTractionN_;

        //- Current patch normal delta
        scalarField deltaN_;

        //- Old time patch normal delta
        scalarField oldDeltaN_;

        //- Effective patch delta, where only positive deltaN is considered
        scalarField deltaEff_;

        //- Unloading patch effective delta: this is the maximum delta effective
        //  that has been reached by the face, so we can monitor if a face
        //  starts to unload
        scalarField unloadingDeltaEff_;

        //- Current energy dissipated in mode-I
        scalarField GI_;

        //- Old time energy dissipated in mode-I
        scalarField oldGI_;

        //- Penalty factor scale
        const scalar penaltyScale_;

        //- Penalty factor for penalty contact
        mutable scalar penaltyFactor_;


    // Private Member functions

        // Calculate normal traction in the damage zone
        vector damageTractionN
        (
            const scalar faceDeltaN,
            const scalar faceSigmaMax,
            const vector& faceN
        ) const;

        //- Calculate and set the penalty factor
        void calcPenaltyFactor() const;

        //- Disallow default bitwise assignment
        void operator=(const modeICohesiveZoneModel&);


public:

    //- Runtime type information
    TypeName("modeI");

    // Static data members


    // Constructors

        //- Construct from dictionary
        modeICohesiveZoneModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict
        );

        //- Construct as copy
        modeICohesiveZoneModel(const modeICohesiveZoneModel&);

        //- Construct and return a clone
        virtual autoPtr<cohesiveZoneModel> clone() const
        {
            return
                autoPtr<cohesiveZoneModel>
                (
                    new modeICohesiveZoneModel(*this)
                );
        }


    // Destructor

        virtual ~modeICohesiveZoneModel();


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

            //- Crack and damage field for visualisation
            //  cracked faces set to 2, damaged but not cracked set to 1
            virtual tmp<scalarField> crackingAndDamage() const;

            //- Write dictionary
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
