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
    variableMixedModeCohesiveZoneModel

Description
    Mixed-mode cohesive zone model where the ratio of normal to shear traction
    components is allowed to varying depending on the ratio of normal to shear
    displacement components.

Author
    Philip Cardif, UCD/UT. All rights reserved.

SourceFiles
    variableMixedModeCohesiveZoneModel.C

\*---------------------------------------------------------------------------*/

#ifndef variableMixedModeCohesiveZoneModel_H
#define variableMixedModeCohesiveZoneModel_H

#include "cohesiveZoneModel.H"
#include "solidCohesiveFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class variableMixedModeCohesiveZoneModel Declaration
\*---------------------------------------------------------------------------*/

class variableMixedModeCohesiveZoneModel
:
    public cohesiveZoneModel
{
    // Private data

        //- Strength in mode-I
        const dimensionedScalar sigmaMax_;

        //- Strength in mode-II
        const dimensionedScalar tauMax_;

        //- Critical mode-I fracture energy
        const dimensionedScalar GIc_;

        //- Critical mode-II fracture energy
        const dimensionedScalar GIIc_;

        //- If a face has cracked:
        //  1 = cracked
        //  0 = damaged but not cracked
        Field<bool> cracked_;

        //- Current normal traction
        scalarField tractionN_;

        //- Old time normal traction
        scalarField oldTractionN_;

        //- Current shear traction
        scalarField tractionS_;

        //- Old time shear traction
        scalarField oldTractionS_;

        //- Current normal delta
        scalarField deltaN_;

        //- Old time normal delta
        scalarField oldDeltaN_;

        //- Current shear delta
        scalarField deltaS_;

        //- Old time shear delta
        scalarField oldDeltaS_;

        //- Effective delta, where only positive deltaN is considered
        scalarField deltaEff_;

        //- Unloading effective delta: this is the maximum delta effective that
        //  has been reached by the face, so we can monitor if a face starts to
        //  unload
        scalarField unloadingDeltaEff_;

        //- Current energy dissipated in mode-I
        scalarField GI_;

        //- Old time energy dissipated in mode-I
        scalarField oldGI_;

        //- Current energy dissipated in mode-II
        scalarField GII_;

        //- Old time energy dissipated in mode-II
        scalarField oldGII_;

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

        // Calculate normal traction in the damage zone
        vector damageTractionN
        (
            const scalar faceDeltaN,
            const scalar faceDeltaS,
            const scalar faceSigmaMax,
            const scalar faceTauMax,
            const vector& faceN
        ) const;

        // Calculate shear traction in the damage zone
        vector damageTractionS
        (
            const scalar faceDeltaN,
            const scalar faceDeltaS,
            const scalar faceSigmaMax,
            const scalar faceTauMax,
            const vector& faceN,
            const vector& faceDelta
        ) const;

        //- Calculate and set the penalty factor
        void calcPenaltyFactor() const;

        //- Disallow default bitwise assignment
        void operator=(const variableMixedModeCohesiveZoneModel&);


public:

    //- Runtime type information
    TypeName("variableMixedMode");

    // Static data members


    // Constructors

        //- Construct from dictionary
        variableMixedModeCohesiveZoneModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict
        );

        //- Construct as copy
        variableMixedModeCohesiveZoneModel
        (
            const variableMixedModeCohesiveZoneModel&
        );

        //- Construct and return a clone
        virtual autoPtr<cohesiveZoneModel> clone() const
        {
            return
                autoPtr<cohesiveZoneModel>
                (
                    new variableMixedModeCohesiveZoneModel(*this)
                );
        }


    // Destructor

        virtual ~variableMixedModeCohesiveZoneModel();


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
