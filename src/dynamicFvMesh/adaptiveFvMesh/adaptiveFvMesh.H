/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 Tyler Voskuilen
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-05-2020  Jeff Heylmun:   Modified original dynamicRefineBalanceFvMesh class
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
    Foam::adaptiveFvMesh

SourceFiles
    adaptiveFvMesh.C
    adaptiveFvMeshTemplates.C

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

#ifndef adaptiveFvMesh_H
#define adaptiveFvMesh_H

#include "dynamicFvMesh.H"
#include "errorEstimator.H"
#include "hexRef.H"
#include "PackedBoolList.H"
#include "Switch.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class adaptiveFvMesh Declaration
\*---------------------------------------------------------------------------*/

class adaptiveFvMesh
:
    public dynamicFvMesh
{
protected:

        autoPtr<errorEstimator> error_;

        //- Mesh cutting engine
        autoPtr<hexRef> meshCutter_;

        //- Dump cellLevel for postprocessing
        Switch dumpLevel_;

        //- Fluxes to map
        HashTable<word> correctFluxes_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;

        //- Number of protected cells
        label nProtected_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCell_;

        //- Does the mesh get balanced
        bool balance_;

        //- Static dictionary for runTime balancing
        dictionary decompositionDict_;

        //- Decomposition method
        autoPtr<decompositionMethod> decomposer_;


    // Protected Member Functions

        label topParentID(const label p) const;

        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCell_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCell) const;

        //- Read the projection parameters from dictionary
        void readDict();


        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);


        // Selection of cells to un/refine

            //- Calculates approximate value for refinement level so
            //  we don't go above maxCell
            scalar getRefineLevel
            (
                const label maxCells,
                const label maxRefinement,
                const scalar refineLevel,
                const scalarField&
            ) const;

            //- Get per cell max of connected point
            scalarField maxPointField(const scalarField&) const;

            //- Get point max of connected cell
            scalarField maxCellField(const volScalarField&) const;

            scalarField cellToPoint(const scalarField& vFld) const;

            scalarField error
            (
                const scalarField& fld,
                const scalar minLevel,
                const scalar maxLevel
            ) const;

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
                const label maxRefinement,
                const PackedBoolList& candidateCell
            ) const;

            //- Select points that can be unrefined.
            virtual labelList selectUnrefinePoints
            (
                const scalar unrefineLevel,
                const PackedBoolList& markedCell,
                const scalarField& pFld
            ) const;

            //- Extend markedCell with cell-face-cell.
            void extendMarkedCells(PackedBoolList& markedCell) const;

            //- Check all cells have 8 anchor points
            void checkEightAnchorPoints
            (
                PackedBoolList& protectedCell,
                label& nProtected
            ) const;

            //- Set protected cells
            void setProtectedCells();

            //- Map non-flux surface<Type>Fields for new internal faces
            //  (from cell splitting)
            template <class T>
            void mapNewInternalFaces(const labelList& faceMap);

            //- Overload update mesh to include other methods
            virtual void updateMesh(const mapPolyMesh& mpm)
            {
                fvMesh::updateMesh(mpm);
            }



public:

    //- Runtime type information
    TypeName("adaptiveFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit adaptiveFvMesh(const IOobject& io);

        //- Disallow default bitwise copy construction
        adaptiveFvMesh(const adaptiveFvMesh&) = delete;


    //- Destructor
    virtual ~adaptiveFvMesh();


    // Member Functions

        //- Map all fields in time using given map
        virtual void mapFields(const mapPolyMesh& mpm);

        //- Direct access to the refinement engine
        const hexRef& meshCutter() const
        {
            return meshCutter_();
        }

        //- Cells which should not be refined/unrefined
        const PackedBoolList& protectedCell() const
        {
            return protectedCell_;
        }

        //- Cells which should not be refined/unrefined
        PackedBoolList& protectedCell()
        {
            return protectedCell_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Template to update all volField boundaries
        template<class Type> void correctBoundaries();

        //- Balance the mesh
        void balance();




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
        void operator=(const adaptiveFvMesh&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#ifdef NoRepository
    #include "adaptiveFvMeshTemplates.C"
#endif

#endif

// ************************************************************************* //
