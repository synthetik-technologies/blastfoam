/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
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
    Foam::fvMeshRefiner

SourceFiles
    fvMeshRefiner.C

Description
    A general class for the selection of mesh refinement method

\*---------------------------------------------------------------------------*/

#ifndef fvMeshRefiner_H
#define fvMeshRefiner_H

#include "polyMeshRefiner.H"
#include "mapDistributePolyMesh.H"
#include "fvMeshBalance.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshRefiner;
typedef MeshObject
<
    fvMesh,
    UpdateableMeshObject,
    fvMeshRefiner
> FvMeshRefiner;

/*---------------------------------------------------------------------------*\
                     Class fvMeshRefiner Declaration
\*---------------------------------------------------------------------------*/

class fvMeshRefiner
:
    public FvMeshRefiner
{
protected:

        //- Reference to the mesh
        fvMesh& mesh_;

        //- polyMesh refiner
        autoPtr<polyMeshRefiner> refiner_;

        //- Balancer
        fvMeshBalance balancer_;

        //- Number of balance steps done so far
        mutable label nBalanceIterations_;

        //- Balancing interval
        label balanceInterval_;

        //- When balancing starts
        scalar beginBalance_;

        //- When balancing ends
        scalar endBalance_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;


        // Saved fields for mapping

            //- Saved volumes as scalarFields
            mutable DimensionedField<scalar, volMesh>* V0OldPtr_;
            mutable DimensionedField<scalar, volMesh>* V00OldPtr_;


    // Protected Member Functions

        //- Preupdate by checking if refinement or unrefinement can occur
        bool preUpdate();

        //- Check if balancing should happen
        virtual bool canBalance(const bool incr = false) const;


public:

    //- Runtime type information
    TypeName("fvMeshRefiner");


        // Declare run-time constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshRefiner,
            fvMesh,
            (fvMesh& mesh),
            (mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            fvMeshRefiner,
            dictionary,
            (fvMesh& mesh, const dictionary dict, const bool force, const bool read),
            (mesh, dict, force, read)
        );


    // Constructors

        //- Construct from mesh
        explicit fvMeshRefiner(fvMesh& mesh);

        //- Construct from mesh and dictionary
        fvMeshRefiner
        (
            fvMesh& mesh,
            const dictionary& dict,
            const bool force = false,
            const bool read = true
        );

        //- Disallow default bitwise copy construction
        fvMeshRefiner(const fvMeshRefiner&) = delete;


    // Selectors

        //- Select, construct and return the fvMeshRefiner
        static autoPtr<fvMeshRefiner> New(fvMesh&);

        static autoPtr<fvMeshRefiner> New
        (
            fvMesh& mesh,
            const dictionary& dict,
            const bool force = false,
            const bool read = true
        );


    //- Destructor
    virtual ~fvMeshRefiner();


    // Member Functions

        // Access

            //- Is the mesh refining
            inline bool isRefining() const
            {
                return refiner_->isRefining();
            }

            //- Is the mesh unrefining
            inline bool isUnrefining() const
            {
                return refiner_->isUnrefining();
            }

            //- Is the mesh being balanced
            inline bool isBalancing() const
            {
                return refiner_->isBalancing();
            }

            //- Return the fvMeshBalancer
            inline const fvMeshBalance& balancer() const
            {
                return balancer_;
            }

            //- Access the fvMeshBalancer
            inline fvMeshBalance& balancer()
            {
                return balancer_;
            }

            //- Return the polyMeshRefiner
            inline const polyMeshRefiner& refiner() const
            {
                return refiner_();
            }

            //- Access the polyMeshRefiner
            inline polyMeshRefiner& refiner()
            {
                return refiner_();
            }

            //- Return the number of buffer layers
            inline label nRefinementBufferLayers() const
            {
                return refiner_->nRefinementBufferLayers_;
            }

            //- Return the number of buffer layers
            inline label nUnrefinementBufferLayers() const
            {
                return refiner_->nUnrefinementBufferLayers_;
            }

            //- Force refinement/unrefinement
            inline void setForce(const bool force)
            {
                refiner_->setForce(force);
            }

            //- Return the cell level
            inline const labelList& cellLevel() const
            {
                return refiner_->cellLevel();
            }

            //- Return the point level
            inline const labelList& pointLevel() const
            {
                return refiner_->pointLevel();
            }

            //- Return the location mapper
            inline const locationMapper& locMapper() const
            {
                return refiner_->locMapper();
            }

            //- Tell the location mapper to store data
            inline void needLocationMap(const bool needMap = true) const
            {
                locMapper().needMap(needMap);
            }


        //- Read the projection parameters from dictionary
        virtual void readDict(const dictionary& parentDict);

        //- Refine the mesh
        virtual bool refine
        (
            const scalarField& error,
            const labelList& maxCellLevel = labelList(),
            const scalar lowerRefineLevel = sqrt(small),
            const scalar upperRefineLevel = great,
            const scalar unrefineLevel = -sqrt(small)
        );

        //- Balance the mesh
        virtual bool balance();

        //- Overload update mesh to include other methods
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Move points
        virtual bool movePoints()
        {
            return false;
        }

        //- Distribute additional information
        virtual void distribute(const mapDistributePolyMesh& map);


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
        void operator=(const fvMeshRefiner&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
