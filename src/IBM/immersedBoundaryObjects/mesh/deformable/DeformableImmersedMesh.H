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
    Foam::DeformableImmersedMesh

Description


SourceFiles
    DeformableImmersedMesh.H
    DeformableImmersedMesh.C

\*---------------------------------------------------------------------------*/

#ifndef DeformableImmersedMesh_H
#define DeformableImmersedMesh_H

#include "immersedMeshMapper.H"
#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dynamicBlastFvMesh;

/*---------------------------------------------------------------------------*\
                      Class DeformableImmersedMesh Declaration
\*---------------------------------------------------------------------------*/
\
template<class ImmersedType>
class DeformableImmersedMesh
:
    public ImmersedType
{
    // Private data

        //- velocity boundary conditions
        dynamicBlastFvMesh& fvMesh_;

        //- Solid mechanics model
        autoPtr<solidModel> solidPtr_;

        //- Mapper to and from the immersedMesh
        immersedMeshMapper* mapperPtr_;

        //- Patch name
        word patchName_;

        autoPtr<pointField> points0Ptr_;


public:

    // Constructors

        //- Construct from constant and state dictionaries
        DeformableImmersedMesh
        (
            const polyMesh& mesh,
            const dictionary& dict,
            const dictionary& stateDict
        );

    //- Destructor
    ~DeformableImmersedMesh();


    // Member Functions

        //- Does this class handle temperature
        bool solveTemperature() const
        {
            return true;
        }

        //- Return the patch name
        virtual word patchName() const
        {
            return patchName_;
        }

        //- Solve
        virtual void solve();

        //- Return pointer to immersed fvMesh
        const fvMesh& immersedMesh() const;

        //- Return a pointer to the immersedMeshMapper
        virtual const immersedMeshMapper* mapper() const;

        //- Move mesh points
        virtual void movePoints();

        //- Initialize
        virtual void initialize();

        //- Return the current Courant number
        virtual scalar CoNum() const;

        //- Return the max Courant number
        virtual scalar maxCoNum() const;

        //- Write
        virtual void write(Ostream&) const;
        virtual void write(dictionary&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DeformableImmersedMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
