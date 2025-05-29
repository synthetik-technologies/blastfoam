/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021-2022 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::fvMeshDistributors::redistributor

Description
    Dynamic mesh redistribution using the distributor specified in
    decomposeParDict

Usage
    Example of single field based refinement in all cells:
    \verbatim
    distributor
    {
        type            distributor;

        libs            ("libblastFvMeshDistributors.so");

        // How often to redistribute
        redistributionInterval  10;

        // Maximum fractional cell distribution imbalance
        // before rebalancing
        maxImbalance    0.1;
    }
    \endverbatim

SourceFiles
    redistributorFvMeshDistributor.C

\*---------------------------------------------------------------------------*/

#ifndef redistributorFvMeshDistributor_H
#define redistributorFvMeshDistributor_H

#include "fvMeshDistributor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    class decompositionMethod;

namespace fvMeshDistributors
{

/*---------------------------------------------------------------------------*\
                          Class redistributor Declaration
\*---------------------------------------------------------------------------*/

class redistributor
:
    public fvMeshDistributor
{
public:
    // Public static functions

        //- Add a constraint
        static dictionary& constraints(dictionary&);

        //- Preserve a faceZone
        static void preserveFaceZone(const wordRe&, dictionary&);

        //- Preserve a faceSet on a single processor
        static void singleProcessorFaceSet
        (
            const word&,
            const label,
            dictionary&
        );

        //- Preserve a patch
        static void preservePatch(const wordRe&, dictionary&);

        //- Preserve baffles
        static void preserveBaffles(dictionary&);


protected:

    // Protected Member Data

        //- Dictionary used for domain decomposition
        dictionary decompositionDict_;

        //- Cache the decomposer/distributor
        autoPtr<decompositionMethod> decomp_;

        //- Is balancing enabled
        bool balance_;

        //- Force balancing, even if other criteria is not met
        bool force_;

        //- Time-step interval between redistribution calls
        label balanceInterval_;

        //- Time to begin balancing
        scalar beginBalance_;

        //- Time to end balancing
        scalar endBalance_;

        //- Maximum imbalance between the ideal number of cells per processor
        //  and the maximum or minimum as a ratio mag(1 - nCells/idealNcells)
        scalar maxImbalance_;

        //- The time index used for updating
        label timeIndex_;

        //- Current iteration (number of times topoChanges have happened)
        label iter_;


    // Protected Member Functions

        //- Set the decomposer
        decompositionMethod& decomp();

        //- Distribute the mesh and mesh data
        void distribute(const labelList& distribution, const bool dist);


public:

    //- Runtime type information
    TypeName("redistributor");


    // Constructors

        //- Construct from fvMesh
        explicit redistributor(fvMesh& mesh, const bool read = true);

        //- Disallow default bitwise copy construction
        redistributor(const redistributor&) = delete;


    //- Destructor
    virtual ~redistributor();


    // Member Functions

        //- Force the distributor to balance
        void setForce(const bool f)
        {
            force_ = f;
        }

        //- Distribute the mesh
        virtual bool update();

        //- force distribute the mesh, optional call to fvMesh::distribute
        bool forceUpdate(const bool dist);

        //- Update corresponding to the given map
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Update corresponding to the given distribution map
        virtual void distribute(const polyDistributionMap&);

        //- Was the given constraint type found?
        bool constraintFound(const word& type) const;

        //- Set the constraint dictionary
        void addConstraint(const word& name, const dictionary& dict);


    // Reading

        //- Read the projection parameters from dictionary
        void readDict();

        //- Read the projection parameters from dictionary
        void readDict(const dictionary& dict);

    // Writing

        //- Write using given format, version and compression
        virtual bool write(const bool write = true) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const redistributor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvMeshDistributors
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
