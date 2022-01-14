/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020 Synthetik Applied Technologies: |   Modified original
                            dynamicRefineBalanceBlastFvMesh class
                            to be more appilcable to compressible flows.
                            Improved compatibility with snappyHexMesh.
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

Authors
    T.G. Voskuilen ( https://github.com/tgvoskuilen/meshBalancing )
    Daniel Deising <deising@mma.tu-darmstadt.de>
    Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
    All rights reserved.

Description
    A fvMesh with run-time load balancing.

    Updated to OpenFOAM-7 from the github version at:
    https://github.com/ElsevierSoftwareX/SOFTX_2018_143.git

    Reference:
    Rettenmaier, D., Deising, D., Ouedraogo, Y., Gjonaj, E., De Gersem, H., Bothe, D., Tropea, C., Marschall, H.
    "Load balanced 2D and 3D adaptive mesh refinement in OpenFOAM" (2019),
    SoftwareX, volume 10, 2352-7110,

    This code has been developed by :
        Daniel Rettenmaier (main developer).

    Method Development and Intellectual Property :
        T.G. Voskuilen (Purdue University)
        Timothée Pourpoint <timothee@purdue.edu> (Purdue University)
        Daniel Rettenmaier <rettenmaier@gsc.tu-darmstadt.de>
        Daniel Deising <deising@mma.tu-darmstadt.de>
        Holger Marschall <marschall@csi.tu-darmstadt.de>
        Dieter Bothe <bothe@csi.tu-darmstadt.de>
        Cameron Tropea <ctropea@sla.tu-darmstadt.de>

        School of Aeronautics and Astronautics Purdue University
        Mathematical Modeling and Analysis
        Institute for Fluid Mechanics and Aerodynamics
        Center of Smart Interfaces
        Technische Universitaet Darmstadt

    If you use this software for your scientific work or your publications,
    please don't forget to acknowledge explicitly the use of it.

    Additional modifications not part of the original work include the use of
    error estimators, improved stability with castellated mesh, and fewer
    required user inputs.

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

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Stored data for access outside of the class

            //- Is refinement/unrefinement forced, i.e. for a meshing process,
            //  not while a case is running
            bool force_;

            // Refinement interval
            label refineInterval_;

            //- When unrefinement starts
            scalar beginUnrefine_;

            //- Maximum number of cells
            label maxCells_;

            //- Number of buffer layers
            //  this is different that the standard dynamicRefineFvMesh since
            //  since this number of cells is refined, not just blocked from
            //  unrefinement
            label nBufferLayers_;

            //- Dump cellLevel for postprocessing
            Switch dumpLevel_;

            //- Protected patches
            wordList protectedPatches_;



        //- Has the mesh been balanced this step
        mutable bool balanced_;

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

        //- Extend markedCell with cell-face-cell.
        void extendMarkedCells
        (
            PackedBoolList& markedCell,
            const labelList& maxRefinement,
            const bool top
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

            //- Return the number of buffer layers
            label nBufferLayers() const
            {
                return nBufferLayers_;
            }

            //- Force refinement/unrefinement
            void setForce(const bool force)
            {
                force_ = force;
            }

            //- Return the cell level
            virtual const labelList& cellLevel() const = 0;


        //- Read the projection parameters from dictionary
        virtual void readDict(const dictionary& parentDict);

        //- Refine the mesh
        virtual bool refine
        (
            const scalarField& error,
            const labelList& maxCellLevel = labelList(),
            const scalar lowerRefineLevel = small,
            const scalar upperRefineLevel = great,
            const scalar unrefineLevel = -small
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
            const fvMesh& mesh,
            labelList& cells,
            labelList& maxCellLevel,
            const label level
        );


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
