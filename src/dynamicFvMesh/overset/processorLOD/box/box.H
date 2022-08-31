/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::processorLODs::box

Description
    Creates the parallel distribution map by describing the source and target
    objects using box shapes.

    A single box is created for the source object, which is then split using
    2x2x2 refinement based on the number of remote target objects that overlap.
    The refinement is local between communicating processor pairs, where the
    refinement continues until a threshold number of remote target objects per
    source box is achieved.

\*---------------------------------------------------------------------------*/

#ifndef processorLODs_box
#define processorLODs_box

#include "processorLOD.H"
#include "treeBoundBox.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace processorLODs
{

/*---------------------------------------------------------------------------*\
                             Class Box Declaration
\*---------------------------------------------------------------------------*/

class box
:
    public processorLOD
{
protected:

    // Protected data

        // Flags to indicate what to do with a box

            //- Drop/discard
            static const label DROP;

            //- Refine
            static const label REFINE;

            //- Fixed - do not touch
            static const label FIXED;


    //- Reference to the source points
    const UList<point>& srcPoints_;

    //- Reference to the target points
    const UList<point>& tgtPoints_;

    //- Per processor, the list of src bound boxes
    //  Note: using treeBoundBox to take advantage of subBbox() method
    List<List<treeBoundBox>> boxes_;

    //- Maximum number of refinement iterations
    label nRefineIterMax_;


    // Caching controls

        //- Number of iterations before element indices are cached
        static const label nStartUpIter;

        //- Addressing per proc of new to old bound boxes
        List<labelList> newToOld_;

        //- Indices of elements in the src boxes
        List<List<labelList>> boxSrcElems_;


    // Private Member Functions

        //- Helper function to write the boxes in OBJ format
        void writeBoxes
        (
            const List<DynamicList<treeBoundBox>>& fixedBoxes,
            const label iter
        ) const;


        virtual boundBox calcSrcBox(const label srcObji) const = 0;
        virtual boundBox calcTgtBox(const label tgtObji) const = 0;

        //- Set the box refinement flags
        void setRefineFlags
        (
            const label refineIter,
            const label nTgtObjects,
            List<labelHashSet>& fixedSendElems,
            List<List<labelList>>& localTgtElems,
            List<labelList>& refineFlags,
            labelList& nElems
        ) const;

        void refineBox
        (
            const label boxi,
            const label refineIter,
            const label nSrcElem,
            const treeBoundBox& origBox,
            DynamicList<treeBoundBox>& procBoxes,
            DynamicList<labelList>& procBoxElems,
            DynamicList<label>& procNewToOld
        ) const;

        void refineBox
        (
            const label boxi,
            const labelList& srcAddr,
            const treeBoundBox& origBox,
            DynamicList<treeBoundBox>& procBoxes,
            DynamicList<labelList>& procBoxElems,
            DynamicList<label>& procNewToOld
        ) const;


        //- Apply the box refinements
        //  \return true if still refining
        bool doRefineBoxes
        (
            const label refineIter,
            const label nSrcFaces,
            const List<labelList>& refineFlags,
            const labelList& nElems,
            List<DynamicList<treeBoundBox>>& fixedBoxes
        );

        autoPtr<mapDistribute> createMap
        (
            const label nSrcElems,
            const label nTgtElems
        );

        //- Use the current list of send elements to create the mapDistribute
        autoPtr<mapDistribute> createLODMap
        (
            List<labelList>& sendElems
        ) const;


public:

    //- Runtime type information
    TypeName("box");

    //- Construct from list of points
    box
    (
        const UList<point>& srcPoints,
        const UList<point>& tgtPoints,
        const label maxObjectsPerLeaf,
        const label nObjectsOfType,
        const label nRefineIterMax = 100
    );

    //- Destructor
    virtual ~box() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace processorLODs

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
