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
    frictionContactModel

Description
    Model for the the contact friction.

    The actual friction law (e.g. Coulomb) is specified by the frictionLaw
    in the frictionLaw.

SourceFiles
    frictionContactModel.C
    newFrictionContactModel.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef frictionContactModel_H
#define frictionContactModel_H

#include "IOdictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "volFields.H"
#include "tmp.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class frictionContactModel Declaration
\*---------------------------------------------------------------------------*/

class frictionContactModel
{
    // Private data

        //- Name
        const word name_;

        //- Reference to primitivePatch
        const fvPatch& patch_;

        //- Master patch index
        const label masterPatchID_;

        //- Slave patch index
        const label slavePatchID_;

        // For visualistion of sticking and slipping faces
        // face not in contact is set to 0.0
        //     slipping face is set to 1.0
        //     sticking face is set to 2.0
        scalarField stickSlipFaces_;


    // Private Member Functions

        //- Disallow copy construct
        //frictionContactModel(const frictionContactModel&);

        //- Disallow default bitwise assignment
        void operator=(const frictionContactModel&);


public:

    //- Runtime type information
    TypeName("frictionContactModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            frictionContactModel,
            dictionary,
            (
                const word name,
                // THIS IS ONLY USED TO ACCESS THE MESH
                // We should instead pass the mesh
                const fvPatch& patch,
                const dictionary& dict,
                const label masterPatchID,
                const label slavePatchID
            ),
            (
                name,
                patch,
                dict,
                masterPatchID,
                slavePatchID
            )
        );


    // Selectors

        //- Return a reference to the selected mechanical model
        static autoPtr<frictionContactModel> New
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );


    // Constructors

        //- Construct from dictionary
        frictionContactModel
        (
            const word& name,
            const fvPatch& patch,
            const dictionary& dict,
            const label masterPatchID,
            const label slavePatchID
        );

        //- Construct as copy
        frictionContactModel(const frictionContactModel&);

        //- Construct and return a clone
        virtual autoPtr<frictionContactModel> clone() const = 0;


    // Destructor

        virtual ~frictionContactModel()
        {}


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return reference to mesh
        const fvMesh& mesh() const
        {
            return patch_.boundaryMesh().mesh();
        }

        //- Correct contatc model
        virtual void correct
        (
            const vectorField& slavePressure,
            const vectorField& slaveFaceNormals,
            const scalarField& slavePatchAreaInContact,
            const vectorField& slaveDU,
            const vectorField& masterDUInterpToSlave
        ) = 0;

        //- It is the start of a new time step
        virtual void newTimeStep() const
        {}

        //- Return slave friction traction
        virtual const vectorField& slaveTraction() const = 0;

        //- Return slave friction traction for master.
        //  This has the same size as the slave patch
        virtual const vectorField& slaveTractionForMaster() const
        {
            return slaveTraction();
        }

        //- Return master patch ID
        virtual label masterPatchID() const
        {
            return masterPatchID_;
        }

        //- Return master patch ID
        virtual label slavePatchID() const
        {
            return slavePatchID_;
        }

        //- Return stick slip faces field
        // virtual volScalarField& stickSlipFaces()
        virtual scalarField& stickSlipFaces()
        {
            return stickSlipFaces_;
        }

        //- Return slip field
        virtual const vectorField& slip() const = 0;

        //- Perform any reauired field mapping e.g. after a topo change
        virtual void autoMap(const fvPatchFieldMapper& m);

        //- Write model dictionary
        virtual void writeDict(Ostream& os) const
        {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
