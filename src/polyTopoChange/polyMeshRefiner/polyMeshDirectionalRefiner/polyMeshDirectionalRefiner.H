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
    Foam::polyMeshDirectionalRefiner

SourceFiles
    polyMeshDirectionalRefiner.C

Description
    Mesh refiner that refines in one direction (x, y, or z) using the
    undoableMeshCutter mesh refiner. This is a simplified implementation of
    the refinemMesh utility. Currently there is no ability to unrefine the
    mesh

\*---------------------------------------------------------------------------*/

#ifndef polyMeshDirectionalRefiner_H
#define polyMeshDirectionalRefiner_H

#include "polyMeshRefiner.H"
#include "newCellCuts.H"
#include "refineCell.H"
#include "undoableMeshCutter.H"
#include "hexCellLooper.H"
#include "geomCellLooper.H"
#include "locationMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class polyMeshDirectionalRefiner Declaration
\*---------------------------------------------------------------------------*/

class polyMeshDirectionalRefiner
:
    public polyMeshRefiner
{
public:
    // Public static functions

        //- Given map from original to added cell set the refineCell for
        //  the added cells to be equal to the one on the original cells.
        static void addCells(const Map<label>&, List<refineCell>&);

        //- Given map from original to added cell set the vectorField for
        //  the added cells to be equal to the one on the original cells.
        static void update(const Map<label>&, vectorField&);

        //- Given map from original to added cell add the added cell to the
        //  list of labels
        static void addCells(const Map<label>&, labelList& labels);

        //- Return a direction given a word
        static vector getDirection(const word& d);


private:
    // Private data

        //- Cell level
        labelIOList cellLevel_;

        //- Point level
        labelIOList pointLevel_;

        //- Directions to refine
        vector refinementDirection_;

        //- Cell walker
        autoPtr<cellLooper> cellWalker_;

        // Construct undoable refinement topology modifier.
        // Note: undoability switched off.
        // Might want to reconsider if needs to be possible. But then can always
        // use other constructor.
        autoPtr<undoableMeshCutter> cutter_;


    // Private Member Functions

        //- Add new cells from map to overall list (addedCells_).
        void addCells
        (
            const Map<label>&,
            labelListList&
        );

        //- Refine cells in cellLabels_ in directions mentioned.
        labelListList refineDirections
        (
            List<vectorField>& cellDirections,
            List<labelList>& cellsToRefine
        );

        //- Refine the given cells
        labelListList refine(const labelList& origCellsToRefine);


public:

    //- Runtime type information
    TypeName("directionalRefiner");


    // Constructors

        //- Construct from polyMesh
        explicit polyMeshDirectionalRefiner(polyMesh& mesh);

        //- Construct from polyMesh and dictionary
        polyMeshDirectionalRefiner
        (
            polyMesh& mesh,
            const dictionary& dict,
            const bool read = true
        );

        //- Disallow default bitwise copy construction
        polyMeshDirectionalRefiner(const polyMeshDirectionalRefiner&) = delete;


    //- Destructor
    virtual ~polyMeshDirectionalRefiner();


    // Member Functions

        // Access

            //- Return the cell level
            virtual const labelList& cellLevel() const
            {
                return cellLevel_;
            }

            //- Return the point level
            virtual const labelList& pointLevel() const
            {
                return pointLevel_;
            }


        //- Read the projection parameters from dictionary
        void readDict(const dictionary& parentDict);

        //- Refine the mesh
        virtual bool refine
        (
            const scalarField& error,
            const labelList& maxCellLevel = labelList(),
            const scalar lowerRefineLevel = small,
            const scalar upperRefineLevel = great,
            const scalar unrefineLevel = -small
        );

        //- Overload update mesh to include other methods
        virtual void topoChange(const polyTopoChangeMap&);

        //- Update from another mesh using the given map
        virtual void mapMesh(const polyMeshMap&);

        //- Distribute additional information
        virtual void distribute(const polyDistributionMap& map);


    // Writing

        //- Write using given format, version and compression
        virtual bool write(const bool write) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const polyMeshDirectionalRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
