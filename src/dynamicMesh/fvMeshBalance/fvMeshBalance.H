/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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
protected:

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Dictionary for runTime balancing
        dictionary decompositionDict_;

        //- Decomposition method
        autoPtr<decompositionMethod> decomposer_;

        fvMeshDistribute distributor_;

        //- Does the mesh get balanced
        bool balance_;

        //- Allowable imbalance
        scalar allowableImbalance_;


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

        //- Template to update all volField boundaries
        template<class GeoField> void correctBoundaries();

        //- Is balancing enabled
        bool balance() const
        {
            return balance_;
        }

        //- Does the mesh need to be balanced
        bool canBalance() const;

        //- Balance the mesh
        autoPtr<mapDistributePolyMesh> distribute();

        //- Helper: push master point data to collocated points
        template<class Type>
        static void pushUntransformedData(const polyMesh&, Field<Type>&);

        //- Read the projection parameters from dictionary
        void read(const dictionary& dict);

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshBalance&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "fvMeshBalanceTemplates.C"
#endif

#endif

// ************************************************************************* //
