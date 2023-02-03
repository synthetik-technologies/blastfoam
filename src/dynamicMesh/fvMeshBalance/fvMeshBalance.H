/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::fvMeshBalance

SourceFiles
    fvMeshBalance.C
    fvMeshBalanceTemplates.C

Description
    Class used to balance a fvMesh

\*---------------------------------------------------------------------------*/

#ifndef fvMeshBalance_H
#define fvMeshBalance_H

#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class decompositionMethod;
/*---------------------------------------------------------------------------*\
                     Class fvMeshBalance Declaration
\*---------------------------------------------------------------------------*/

class fvMeshBalance
{
public:

    // Global flag if the mesh is being balanced
    static bool balancing;


protected:

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Dictionary for runTime balancing
        mutable dictionary decompositionDict_;

        //- Has the decomposition dict been modified
        mutable bool modified_;

        //- Constraints dictionary
        dictionary* constraintsDict_;

        //- Preserve faceZones dict
        dictionary* preserveFaceZonesDict_;

        //- Preserve faceSets dict
        dictionary* singleProcessorFaceSetsDict_;

        //- Preserve patches dict
        dictionary* preservePatchesDict_;

        //- Preserve baffles dict
        dictionary* preserveBafflesDict_;

        //- Decomposition method
        mutable autoPtr<decompositionMethod> decomposer_;

        fvMeshDistribute distributor_;

        //- Does the mesh get balanced
        bool balance_;

        //- Allowable imbalance
        scalar allowableImbalance_;

        //- Distribution
        mutable labelList distribution_;


    // Private member functions

        //- Set the decomposer
        void makeDecomposer() const;

        //- Check in an internal patch is found
        void checkForInternal() const;


public:

    //- Runtime type information
    TypeName("fvMeshBalance");


    // Constructors

        //- Construct from mesh
        explicit fvMeshBalance(fvMesh& mesh);

        //- Construct from IOobject
        fvMeshBalance(fvMesh& mesh, const dictionary& dict);

        //- Disallow default bitwise copy construction
        fvMeshBalance(const fvMeshBalance&) = delete;


    //- Destructor
    virtual ~fvMeshBalance();


    // Member Functions

        //- Is balancing enabled
        bool balance() const
        {
            return balance_;
        }

        //- Access if balancing enabled
        bool& balance()
        {
            return balance_;
        }

        //- Does the mesh need to be balanced
        bool canBalance() const;

        // Decomposition constraints

            //- Add a constraint
            void addConstraint(const word& dictName, const dictionary& dict);

            //- Preserve a faceZone
            void preserveFaceZone(const wordRe& zoneName);

            //- Preserve a faceSet on a single processor
            void singleProcessorFaceSet
            (
                const word& setName,
                const label proc
            );

            //- Preserve a patch
            void preservePatch(const wordRe& patchName);

            //- Preserve baffles
            void preserveBaffles();

        //- Access the decompositionMethod
        decompositionMethod& decomposer() const;

        //- Balance the mesh
        autoPtr<mapDistributePolyMesh> distribute();

        //- Read the projection parameters from dictionary
        void read(const dictionary& dict);

        //- Write the decompositionDict
        bool write(const bool write = true) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshBalance&) = delete;


    // Static Member Functions

        //- Template to update all processor boundaries
        template<class GeoField>
        static void correctProcessorBoundaries(const fvMesh&);

        //- Helper: push master point data to collocated points
        template<class Type>
        static void pushUntransformedData(const polyMesh&, Field<Type>&);

        // Return is the mesh is being balanced
        static bool isBalancing();
        static bool isBalancing(const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "fvMeshBalanceTemplates.C"
#endif

#endif

// ************************************************************************* //
