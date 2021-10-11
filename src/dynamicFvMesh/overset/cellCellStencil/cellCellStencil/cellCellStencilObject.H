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
    Foam::cellCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef cellCellStencilObject_H
#define cellCellStencilObject_H

#include "cellCellStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cellCellStencilObject;
typedef MeshObject
<
    fvMesh,
    Foam::MoveableMeshObject,
    cellCellStencilObject
> Stencil;

/*---------------------------------------------------------------------------*\
                  Class cellCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class cellCellStencilObject
:
    public Stencil,
    public cellCellStencil
{
    // Private data

        autoPtr<cellCellStencil> stencilPtr_;


public:

    TypeName("cellCellStencilObject");

    // Constructors

        //- Construct with mesh
        explicit cellCellStencilObject
        (
            const fvMesh& mesh,
            const bool update = true
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::MoveableMeshObject,
                cellCellStencilObject
            >(mesh),
            cellCellStencil(mesh),
            stencilPtr_
            (
                cellCellStencil::New
                (
                    mesh,
                    mesh.schemesDict().subDict
                    (
                        "oversetInterpolation"
                    ),
                    update
                )
            )
        {}


    //- Destructor
    virtual ~cellCellStencilObject() = default;


    // Member Functions

        //- Callback for geometry motion
        virtual bool movePoints()
        {
            return stencilPtr_().update();
        }

        //- Update stencils. Return false if nothing changed.
        virtual bool update()
        {
            return stencilPtr_().update();
        }

        //- Return the cell type list
        virtual const labelUList& cellTypes() const
        {
            return stencilPtr_().cellTypes();
        }

        //- Indices of interpolated cells
        virtual const labelUList& interpolationCells() const
        {
            return stencilPtr_().interpolationCells();
        }

        //- Return a communication schedule
        virtual const mapDistribute& cellInterpolationMap() const
        {
            return stencilPtr_().cellInterpolationMap();
        }

        //- Per interpolated cell the neighbour cells (in terms of slots as
        //  constructed by above cellInterpolationMap) to interpolate
        virtual const labelListList& cellStencil() const
        {
            return stencilPtr_().cellStencil();
        }

        //- Weights for cellStencil
        virtual const List<scalarList>& cellInterpolationWeights() const
        {
            return stencilPtr_().cellInterpolationWeights();
        }

        //- Per interpolated cell the interpolation factor. (0 = use
        //  calculated, 1 = use interpolated)
        virtual const scalarList& cellInterpolationWeight() const
        {
            return stencilPtr_().cellInterpolationWeight();
        }

        //- Calculate weights for a single acceptor
        virtual void stencilWeights
        (
            const point& sample,
            const pointList& donorCcs,
            scalarList& weights
        ) const
        {
            stencilPtr_().stencilWeights(sample, donorCcs, weights);
        }

        //- Return the names of any (stencil or mesh specific) fields that
        //  should not be interpolated
        virtual const wordHashSet& nonInterpolatedFields() const
        {
            return stencilPtr_().nonInterpolatedFields();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
