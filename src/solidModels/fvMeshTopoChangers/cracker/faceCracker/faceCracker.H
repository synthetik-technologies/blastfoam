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
    faceCracker

Description
    Face cracker mesh modifier. This modifier takes a set of
    internal face labels and converts them into boundary faces.

SourceFiles
    faceCracker.C
    detachFaceCracker.C

\*---------------------------------------------------------------------------*/

#ifndef faceCracker_H
#define faceCracker_H

#include "polyMesh.H"
#include "polyTopoChange.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// PC, Jun-2015
typedef HashTable
<
    nil,
    labelPair,
    labelPair::Hash<>
> labelPairHashSet;

/*---------------------------------------------------------------------------*\
                          Class faceCracker Declaration
\*---------------------------------------------------------------------------*/

class faceCracker
{
    // Private data

        //- Reference to polyMesh
        const polyMesh& mesh_;

        //- Crack face zone name
        word crackZone_;

        //- Crack patch Name
        word crackPatch_;

        //- List of boundary faces to move to the crack patch
        mutable labelList coupledFacesToBreak_;

        //- Topology change trigger
        mutable bool trigger_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceCracker(const faceCracker&);

        //- Disallow default bitwise assignment
        void operator=(const faceCracker&);

        //- Detach internal faces
        void detachInternalFaces(polyTopoChange&) const;

        //- Detach coupled faces
        void detachCoupledFaces(polyTopoChange&) const;

public:

    //- Runtime type information
    TypeName("faceCracker");


    // Constructors

        //- Construct from components
        faceCracker
        (
            const polyMesh& mesh,
            const word& faceZoneName,
            const word& crackPatchName
        );

        //- Construct from dictionary
        faceCracker
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor

        virtual ~faceCracker();


    // Member Functions

        //- Return master patch ID
        const word& crackPatch() const
        {
            return crackPatch_;
        }

        //- Set faces to break, where an internal face is converted into two
        //  boundary faces in the crack patch, or two coupled faces are moved
        //  to the crack patch
        void setBreak
        (
            polyMesh& mesh,
            const labelList& facesToBreak,
            const boolList& faceFlip,
            const labelList& coupledFacesToBreak
        );

        //- Check for topology change
        virtual bool changeTopology() const;

        //- Set topological change
        virtual void setRefinement(polyTopoChange&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
