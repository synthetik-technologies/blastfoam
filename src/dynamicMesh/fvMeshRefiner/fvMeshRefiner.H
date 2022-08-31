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

#include "hexRef.H"
#include "PackedBoolList.H"
#include "Switch.H"
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

        //- Dictionary used to read parameters
        dictionary dict_;

        //- Balancer
        fvMeshBalance balancer_;

        // Track number of times refinement/unrefinement/balancing is done

            //- Number of refinement steps done so far.
            mutable label nRefinementIterations_;

            //- Number of unrefinement steps done so far.
            mutable label nUnrefinementIterations_;

            //- Number of balance steps done so far
            mutable label nBalanceIterations_;


        //- Is refinement/unrefinement/balancing forced, i.e. for a meshing
        // process
        bool force_;

        //- Does refinement occur
        bool refine_;

        //- Does unrefinement occur
        bool unrefine_;

        //- Does balancing occur
        // Handled in fvMeshBalance
        // bool balance_;


        //- Refinement interval
        label refineInterval_;

        //- Unrefinement interval
        label unrefineInterval_;

        //- Balancing interval
        label balanceInterval_;


        //- When refinement starts
        scalar beginRefine_;

        //- When unrefinement starts
        scalar beginUnrefine_;

        //- When balancing starts
        scalar beginBalance_;


        //- When refinement ends
        scalar endRefine_;

        //- When unrefinement ends
        scalar endUnrefine_;

        //- When balancing ends
        scalar endBalance_;


        //- Maximum number of cells
        label maxCells_;

        //- Number of buffer layers for refinement
        label nRefinementBufferLayers_;

        //- Number of buffer layers for refinement
        label nUnrefinementBufferLayers_;

        //- Protected patches
        wordList protectedPatches_;


        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;




        // Track the state of the refiner

            //- Is the mesh in the process of being refined
            mutable bool isRefining_;

            //- Is the mesh in the process of being unrefined
            mutable bool isUnrefining_;

            //- Is the mesh in the process of being balanced
            mutable bool isBalancing_;


        // Saved fields for mapping

            //- Saved volumes as scalarFields
            mutable DimensionedField<scalar, volMesh>* V0OldPtr_;
            mutable DimensionedField<scalar, volMesh>* V00OldPtr_;


    // Protected Member Functions

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Select candidate cells for refinement
        virtual void selectRefineCandidates
        (
            const scalar lowerRefineLevel,
            const scalar upperRefineLevel,
            const scalarField& vFld,
            PackedBoolList& candidateCell
        ) const;

        //- Subset candidate cells for refinement
        virtual labelList selectRefineCells
        (
            const label maxCells,
            const labelList& maxRefinement,
            const PackedBoolList& candidateCell
        ) const;

        //- Get per cell max of connected point
        scalarField maxPointField(const scalarField&) const;

        //- Get point max of connected cell
        scalarField maxCellField(const scalarField&) const;

        scalarField cellToPoint(const scalarField& vFld) const;

        scalarField error
        (
            const scalarField& fld,
            const scalar minLevel,
            const scalar maxLevel
        ) const;

        //- Preupdate by checking if refinement or unrefinement can occur
        bool preUpdate();

        //- Set the maximum cell level
        void setMaxCellLevel(labelList& maxCellLevel) const;

        //- Check if refinement should happen
        virtual bool canRefine(const bool incr = false) const;

        //- Check if unrefinement should happen
        virtual bool canUnrefine(const bool incr = false) const;

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
            bool isRefining() const
            {
                return isRefining_;
            }

            //- Is the mesh unrefining
            bool isUnrefining() const
            {
                return isUnrefining_;
            }

            //- Is the mesh being balanced
            bool isBalancing() const
            {
                return isBalancing_;
            }

            //- Return the number of buffer layers
            label nRefinementBufferLayers() const
            {
                return nRefinementBufferLayers_;
            }

            //- Return the number of buffer layers
            label nUnrefinementBufferLayers() const
            {
                return nUnrefinementBufferLayers_;
            }

            //- Force refinement/unrefinement
            void setForce(const bool force)
            {
                force_ = force;
            }

            //- Return the cell level
            virtual const labelList& cellLevel() const = 0;

            //- Return the point level
            virtual const labelList& pointLevel() const = 0;

            //- Return the location mapper
            virtual const locationMapper& locMapper() const = 0;

            //- Tell the location mapper to store data
            void needLocationMap(const bool needMap = true) const
            {
                locMapper().needMap(needMap);
            }

            //- Return the fvMeshBalancer
            const fvMeshBalance& balancer() const
            {
                return balancer_;
            }

            //- Access the fvMeshBalancer
            fvMeshBalance& balancer()
            {
                return balancer_;
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
        ) = 0;

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

        static void extendMaxCellLevel
        (
            const polyMesh& mesh,
            labelList& cells,
            labelList& maxCellLevel,
            const label level
        );

        void extendMarkedCells
        (
            PackedBoolList& markedCells,
            const labelList& maxCellLevel,
            const bool isTop,
            const bool force
        );

        //- Extend marked cells across faces given a bool list of already marked
        //  cells
        void extendMarkedCellsAcrossFaces
        (
            PackedBoolList& markedCells
        );

        //- Extend marked cells across points given a bool list of already
        //  marked cells
        void extendMarkedCellsAcrossPoints
        (
            PackedBoolList& markedCells
        );


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


// Helper class for accessing max cell level of faces accross processor patches
template<class Type>
class combineMaxOp
{
    public:
    void operator()(Type& x, const Type& y) const
    {
        x = max(x, y);
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
