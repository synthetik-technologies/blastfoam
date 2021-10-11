/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::cellCellStencil

Description
    Calculation of interpolation stencils.

    Looks up zoneID labelIOList to give the zoning. Wrapped in
    MeshObject as cellCellStencilObject. Kept separate so meshes can
    implement more clever methods (e.g. solid body motion does not require
    full recalculation)

SourceFiles
    cellCellStencil.C
    cellCellStencilObject.C

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencil_H
#define cellCellStencil_H

#include "scalarList.H"
#include "mapDistribute.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapDistribute;

typedef UList<point> pointUList;
typedef List<point> pointList;
typedef SubList<point> pointSubList;
typedef List<pointList> pointListList;


/*---------------------------------------------------------------------------*\
                     Class cellCellStencil Declaration
\*---------------------------------------------------------------------------*/

class cellCellStencil
{
public:

    enum patchCellType
    {
        OTHER = 0,          // not on special patch
        PATCH = 1,          // next to (non-coupled) boundary
        OVERSET = 2         // next to 'overset' boundary
    };

    enum cellType
    {
        CALCULATED = 0,     // normal operation
        INTERPOLATED = 1,   // interpolated
        HOLE = 2            // hole
    };
    static const NamedEnum<cellType, 3> names;


protected:

    // Protected data

        //- Mode type names
        static const NamedEnum<cellType, 3> cellTypeNames_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Set of fields that should not be interpolated
        wordHashSet nonInterpolatedFields_;


    // Protected Member Functions

        //- Count occurrences (in parallel)
        static labelList count(const label size, const labelUList& lst);

        //- Helper: create volScalarField for postprocessing.
        template<class Type>
        static tmp<volScalarField> createField
        (
            const fvMesh& mesh,
            const word& name,
            const UList<Type>&
        );


private:

    // Private Member Functions

        //- No copy construct
        cellCellStencil(const cellCellStencil&) = delete;

        //- No copy assignment
        void operator=(const cellCellStencil&) = delete;


public:

    //- Runtime type information
    TypeName("cellCellStencil");


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellCellStencil,
            mesh,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const bool update
            ),
            (mesh, dict, update)
        );


    // Constructors

        //- Construct from fvMesh
        cellCellStencil(const fvMesh&);

        //- New function which constructs and returns pointer to a
        //  cellCellStencil
        static autoPtr<cellCellStencil> New
        (
            const fvMesh&,
            const dictionary& dict,
            const bool update = true
        );


    //- Destructor
    virtual ~cellCellStencil();


    // Member Functions

        //- Update stencils. Return false if nothing changed.
        virtual bool update() = 0;

        //- Return the cell type list
        virtual const labelUList& cellTypes() const = 0;

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const = 0;

        //- Return a communication schedule
        virtual const mapDistribute& cellInterpolationMap() const = 0;

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const = 0;

        //- Weights for cellStencil
        virtual const List<scalarList>& cellInterpolationWeights() const = 0;

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const = 0;

        //- Calculate weights for a single acceptor
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const = 0;

        //- Return the names of any (stencil or mesh specific) fields that
        //  should not be interpolated
        virtual const wordHashSet& nonInterpolatedFields() const;

        //- Return non-const non-interpolating fields
        virtual wordHashSet& nonInterpolatedFields();

        //- Helper: is stencil fully local
        bool localStencil(const labelUList&) const;

        //- Helper: get reference to registered zoneID. Loads volScalarField
        //  if not registered.
        static const labelIOList& zoneID(const fvMesh&);

        //- Helper: get reference to registered zoneID. Loads volScalarField
        //  if not registered.
        const labelIOList& zoneID() const
        {
            return zoneID(mesh_);
        }

        //- Helper: create cell-cell addressing in global numbering
        static void globalCellCells
        (
            const globalIndex& gi,
            const polyMesh& mesh,
            const boolList& isValidDonor,
            const labelList& selectedCells,
            labelListList& cellCells,
            pointListList& cellCellCentres
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "cellCellStencilTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
