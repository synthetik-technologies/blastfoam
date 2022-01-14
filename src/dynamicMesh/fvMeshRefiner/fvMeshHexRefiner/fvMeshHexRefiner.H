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
    Foam::fvMeshHexRefiner

SourceFiles
    fvMeshHexRefiner.C

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

#ifndef fvMeshHexRefiner_H
#define fvMeshHexRefiner_H

#include "fvMeshRefiner.H"
#include "hexRef.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshHexRefiner Declaration
\*---------------------------------------------------------------------------*/

class fvMeshHexRefiner
:
    public fvMeshRefiner
{
protected:

        //- Mesh cutting engine
        autoPtr<hexRef> meshCutter_;

        //- Number of protected cells
        label nProtected_;

        //- Protected cells (usually since not hexes)
        PackedBoolList protectedCell_;

        //- Is the mesh in the process of being refined
        mutable bool isRefining_;

        //- Is the mesh in the process of being unrefined
        mutable bool isUnrefining_;


    // Protected Member Functions

        //- Calculate cells that cannot be refined since would trigger
        //  refinement of protectedCell_ (since 2:1 refinement cascade)
        void calculateProtectedCells(PackedBoolList& unrefineableCell) const;


        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);

        //- Subset candidate cells for refinement
        virtual labelList selectRefineCells
        (
            const label maxCells,
            const labelList& maxRefinement,
            const PackedBoolList& candidateCell
        ) const;

        //- Check all cells have 8 anchor points
        void checkEightAnchorPoints
        (
            PackedBoolList& protectedCell,
            label& nProtected
        ) const;

        //- Set protected cells
        void setProtectedCells();


public:

    //- Runtime type information
    TypeName("hexRefiner");


    // Constructors

        //- Construct from mesh
        explicit fvMeshHexRefiner
        (
            fvMesh& mesh
        );

        //- Construct from mesh and dictionary
        explicit fvMeshHexRefiner
        (
            fvMesh& mesh,
            const dictionary& dict,
            const bool force = false,
            const bool read = true
        );

        //- Disallow default bitwise copy construction
        fvMeshHexRefiner(const fvMeshHexRefiner&) = delete;


    //- Destructor
    virtual ~fvMeshHexRefiner();


    // Member Functions

        //- Direct access to the refinement engine
        const hexRef& meshCutter() const
        {
            return meshCutter_();
        }

        //- Direct access to the refinement engine
        hexRef& meshCutter()
        {
            return meshCutter_();
        }

        //- Return the cell level
        virtual const labelList& cellLevel() const
        {
            return meshCutter_->cellLevel();
        }

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
        virtual void updateMesh(const mapPolyMesh& mpm);

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
        void operator=(const fvMeshHexRefiner&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
