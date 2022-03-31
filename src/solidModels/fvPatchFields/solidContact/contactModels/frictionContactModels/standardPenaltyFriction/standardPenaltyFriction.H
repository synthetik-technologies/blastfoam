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
    standardPenaltyFriction

Description
    Friction contact law using the penalty method to enforce the friction
    constraints.

SourceFiles
    standardPenaltyFriction.C

Author
    Philip Cardiff, UCD. All rights reserved.
    Matthias Schnaithmann, Universitat Stuttgart. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef standardPenaltyFriction_H
#define standardPenaltyFriction_H

#include "frictionContactModel.H"
#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class standardPenaltyFriction Declaration
\*---------------------------------------------------------------------------*/

class standardPenaltyFriction
:
    public frictionContactModel
{
    // Private data

         //- The contact model dictionary
        dictionary frictionContactModelDict_;

         //- Friction law
        autoPtr<frictionLaw> frictionLawPtr_;

         //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Slave traction
        volVectorField slaveTractionVolField_;

        //- Previous slave traction used for under-relaxation
        //vectorField prevSlaveTraction_;

        //- Slip
        vectorField slip_;

        //- Frictional penalty factor and scaling factor
        scalar frictionPenaltyFactor_;

        //- Frictional penalty scaling factor
        scalar frictionPenaltyScale_;

        //- Under-relaxation factor
        scalar relaxFac_;

        //- Contact iteration number
        label contactIterNum_;


    // Private Member Functions

        //- Return penalty factor
        virtual scalar frictionPenaltyFactor();

        //- Disallow default bitwise copy construct
        //standardPenaltyFriction(const standardPenaltyFriction&);

        //- Disallow default bitwise assignment
        void operator=(const standardPenaltyFriction&);

public:

    //- Runtime type information
    TypeName("standardPenalty");


    // Constructors

        //- Construct from dictionary
        standardPenaltyFriction
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );

        //- Construct as copy
        standardPenaltyFriction(const standardPenaltyFriction&);

        //- Construct and return a clone
        virtual autoPtr<frictionContactModel> clone() const
        {
            return
                autoPtr<frictionContactModel>
                (
                    new standardPenaltyFriction(*this)
                );
        }


    //- Destructor
    virtual ~standardPenaltyFriction();


    // Member Functions

        //- Correct the friction model
        virtual void correct
        (
            const vectorField& slavePressure,
            const vectorField& slaveFaceNormals,
            const scalarField& slavePatchAreaInContact,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        );

        //- Return const reference to slave friction traction
        virtual const vectorField& slaveTraction() const
        {
            return slaveTractionVolField_.boundaryField()[slavePatchID()];
        }

        //- Return non-const reference to slave friction traction
        virtual vectorField& slaveTraction()
        {
            return slaveTractionVolField_.boundaryFieldRef()[slavePatchID()];
        }

        //- calculate penalty factor
        void calcFrictionPenaltyFactor();

        //- Return slip field
        virtual const vectorField& slip() const
        {
            // Note: slip is updated inside the correct function
            return slip_;
        }

        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
