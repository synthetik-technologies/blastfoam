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
    cohesiveZoneInitiation

Description
    Failure class that selects internal mesh faces to break based on the
    initiation criteria of the current cohesive zone model.
    The cohesive zone model is looked up and its 'initiationTractionFraction'
    function is called to determine which faces to break.

    In essense, the current law is a wrapper class that asked the
    cohesiveZoneModel class what faces it would like to break.

Author
    Philip Cardif, UCD/UT. All rights reserved.

SourceFiles
    cohesiveZoneInitiation.C

\*---------------------------------------------------------------------------*/

#ifndef cohesiveZoneInitiation_H
#define cohesiveZoneInitiation_H

#include "faceBreakerLaw.H"
#include "crackPathLimiter.H"
#include "cohesiveZoneModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class cohesiveZoneInitiation Declaration
\*---------------------------------------------------------------------------*/

class cohesiveZoneInitiation
:
    public faceBreakerLaw
{
    // Private data

        // Patch index of the cohesive patch
        mutable label* cohesivePatchIDPtr_;

        // Allow cracks along coupled faces
        const Switch allowCoupledFaces_;

        // List of internal faces to break
        mutable labelList* facesToBreakPtr_;

        // List of coupled faces to break
        mutable labelList* coupledFacesToBreakPtr_;

        // List of tractions of faces to break
        mutable List<vector>* facesToBreakTractionsPtr_;

        // List of tractions of coupled faces to break
        mutable List<vector>* coupledFacesToBreakTractionsPtr_;

        // List of normals of faces to break
        mutable List<vector>* facesToBreakNormalsPtr_;

        // List of normals of coupled faces to break
        mutable List<vector>* coupledFacesToBreakNormalsPtr_;

        //- Crack path limiter to limit the cracks to certain regions
        mutable autoPtr<crackPathLimiter> pathLimiterPtr_;


    // Private Member Functions

        //- Calculate the cohesive patch ID
        void calcCohesivePatchID() const;

        //- Return cohesive patch ID
        label cohesivePatchID() const;

        //- Lookup and return a const reference to the cohesiveZoneModel in the
        //  solidCohesive patch
        const cohesiveZoneModel& cohesiveZone() const;

        //- Calculate all faces to break
        void calcAllFacesToBreak() const;

        //- Disallow default bitwise copy construct
        cohesiveZoneInitiation(const cohesiveZoneInitiation&);

        //- Disallow default bitwise assignment
        void operator=(const cohesiveZoneInitiation&);


public:

    //- Runtime type information
    TypeName("cohesiveZoneInitiation");

    // Static data members


    // Constructors

        //- Construct from dictionary
        cohesiveZoneInitiation
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~cohesiveZoneInitiation();


    // Member Functions

        // Return list of internal faces to break
        virtual const labelList& facesToBreak() const;

        // Return list of coupled faces to break
        virtual const labelList& coupledFacesToBreak() const;

        // Return list of internal faces to break tractions
        virtual const List<vector>& facesToBreakTractions() const;

        // Return list of coupled faces to break tractions
        virtual const List<vector>& coupledFacesToBreakTractions() const;

        // Return list of internal faces to break normals
        virtual const List<vector>& facesToBreakNormals() const;

        // Return list of coupled faces to break normals
        virtual const List<vector>& coupledFacesToBreakNormals() const;

        //- Delete demand driven data
        virtual void clearOut();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
