/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::NondeformableImmersedMesh

Description


SourceFiles
    NondeformableImmersedMeshI.H
    NondeformableImmersedMesh.C
    NondeformableImmersedMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef NondeformableImmersedMesh_H
#define NondeformableImmersedMesh_H

#include "fvMesh.H"
#include "immersedMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class NondeformableImmersedMesh Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class NondeformableImmersedMesh
:
    public ImmersedType
{
    // Private data

        //- velocity boundary conditions
        fvMesh& fvMesh_;

        //- Mapper to and from the immersedMesh
        immersedMeshMapper* mapperPtr_;

        //- Original mesh points
        pointField points0_;

        //- Patch name
        word patchName_;


public:

    // Constructors

        //- Construct from constant and state dictionaries
        NondeformableImmersedMesh
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~NondeformableImmersedMesh();


    // Member Functions

        //- Does this class handle temperature
        bool solveTemperature() const
        {
            return false;
        }

        //- Return the patch name
        virtual word patchName() const
        {
            return patchName_;
        }

        //- Return pointer to immersed fvMesh
        const fvMesh& immersedMesh() const;

        //- Return a pointer to the immersedMeshMapper
        virtual const immersedMeshMapper* mapper() const;

        //- Move mesh points
        virtual void movePoints();

        //- Initialize
        virtual void initialize();

        //- Write
        virtual void write(Ostream&) const;
        virtual void write(dictionary&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "NondeformableImmersedMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
