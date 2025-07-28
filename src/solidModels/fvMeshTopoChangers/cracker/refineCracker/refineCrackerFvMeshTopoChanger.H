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
    crackerFvMesh

Description
    A mesh allowing topological detachment on internal faces of the mesh, based
    on a run-time selectable faceBreakerLaw.

SourceFiles
    refineCrackerFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef refineCrackerFvMesh_H
#define refineCrackerFvMesh_H

#include "crackerFvMesh.H"
#include "errorEstimator.H"
#include "fvMeshRefiner.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class RefineCrackerFvMesh Declaration
\*---------------------------------------------------------------------------*/

template<class CrackerMesh>
class RefineCrackerFvMesh
:
    public CrackerMesh
{
    // Private data

        //- Error Estimator
        autoPtr<errorEstimator> error_;

        //- Mesh cutting engine
        autoPtr<fvMeshRefiner> refiner_;

        label curTimeIndex_;


        // Protected member functions

            //- Read from dynamicFvMesh
            void readDict();

            //- Distribute additional information
            virtual void distribute(const mapDistributePolyMesh& map);


public:

    //- Runtime type information
    TypeName("CrackerFvMesh")


    // Constructors

        //- Construct from IOobject
        explicit RefineCrackerFvMesh(const IOobject& io);


    // Destructor

        virtual ~RefineCrackerFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const fvMeshRefiner& refiner() const
        {
            return refiner_();
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Refine the mesh
        virtual bool refine();


    // Writing

        //- Write using given format, version and compression
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write = true
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const RefineCrackerFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "refineCrackerFvMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
