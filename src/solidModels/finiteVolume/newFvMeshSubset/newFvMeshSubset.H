/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::newFvMeshSubset

Description
    Post-processing mesh subset tool.  Given the original mesh and the
    list of selected cells, it creates the mesh consisting only of the
    desired cells, with the mapping list for points, faces, and cells.

    Puts all exposed internal faces into either
    - a user supplied patch
    - a newly created patch "oldInternalFaces"

    - setCellSubset is for small subsets. Uses Maps to minimize memory.
    - setLargeCellSubset is for largish subsets (>10% of mesh).
      Uses labelLists instead.

    - setLargeCellSubset does coupled patch subsetting as well. If it detects
      a face on a coupled patch 'losing' its neighbour it will move the
      face into the oldInternalFaces patch.

    - if a user supplied patch is used the mapping becomes a problem.
    Do the new faces get the value of the internal face they came from?
    What if e.g. the user supplied patch is a fixedValue 0? So for now
    they get the face of existing patch face 0.

SourceFiles
    newFvMeshSubset.C

\*---------------------------------------------------------------------------*/

#ifndef newFvMeshSubset_H
#define newFvMeshSubset_H

#include "fvMesh.H"
#include "pointMesh.H"
#include "fvPatchFieldMapper.H"
// #include "pointPatchFieldMapper.H"
#include "GeometricField.H"
#include "HashSet.H"
#include "surfaceMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class newFvMeshSubset Declaration
\*---------------------------------------------------------------------------*/

class newFvMeshSubset
:
    public regIOobject
{

public:

    //- Patch-field subset interpolation class
    class patchFieldSubset
    :
        public fvPatchFieldMapper
    {
        label sizeBeforeMapping_;

        labelField directAddressing_;

    public:

        // Constructors

            //- Construct from components
            patchFieldSubset(const label sbm, const labelList& da)
            :
                sizeBeforeMapping_(sbm),
                directAddressing_(da)
            {}

            //- Construct given addressing
            patchFieldSubset(const newFvMeshSubset& ms, const label patchNo)
            :
                sizeBeforeMapping_
                (
                    ms.baseMesh().boundary()[ms.patchMap()[patchNo]].size()
                ),
                directAddressing_
                (
                    static_cast<const labelField&>
                    (
                        labelField::subField
                        (
                            ms.faceMap(),
                            ms.subMesh().boundary()[patchNo].size(),
                            ms.subMesh().boundary()[patchNo].patch().start()
                        )
                    )
                  - ms.baseMesh().boundary()
                    [ms.patchMap()[patchNo]].patch().start()
                )
            {
                // If patchNo supplied so exposed internal faces or uncoupled
                // patch faces get into existing patch what to do with
                // mapping? Truncate to 0 for now.

                const label sz =
                    ms.baseMesh().boundary()[ms.patchMap()[patchNo]].size();

                forAll(directAddressing_, i)
                {
                    label& address = directAddressing_[i];

                    if (address < 0 || address >= sz)
                    {
                        address = 0;
                    }
                }
            }


        // Destructor

            virtual ~patchFieldSubset()
            {}


        // Member Functions

            label size() const
            {
                return directAddressing_.size();
            }

            label sizeBeforeMapping() const
            {
                return directAddressing_.size();
            }

            bool direct() const
            {
                return true;
            }

            const unallocLabelList& directAddressing() const
            {
                return directAddressing_;
            }
    };


    //- Patch-field subset interpolation class
    class pointPatchFieldSubset
    :
        public PointPatchFieldMapper
    {
        //- Addressing
        const labelList& directAddressing_;


    public:

        // Constructors

            //- Construct given addressing
            pointPatchFieldSubset(const labelList& directAddressing)
            :
                directAddressing_(directAddressing)
            {}

        // Destructor

            virtual ~pointPatchFieldSubset()
            {}


        // Member Functions

            label size() const
            {
                return directAddressing_.size();
            }

            label sizeBeforeMapping() const
            {
                return directAddressing_.size();
            }

            bool direct() const
            {
                return true;
            }

            const unallocLabelList& directAddressing() const
            {
                return directAddressing_;
            }
    };


private:

    // Private data

        //- Reference to original mesh
        const fvMesh& baseMesh_;

        //- Subset fvMesh pointer
        fvMesh* newFvMeshSubsetPtr_;

        //- Subset pointMesh pointer
        mutable pointMesh* pointMeshSubsetPtr_;

        //- Point mapping array
        labelList pointMap_;

        //- Face mapping array
        labelList faceMap_;

        //- Cell mapping array
        labelList cellMap_;

        //- Patch mapping array
        labelList patchMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        newFvMeshSubset(const newFvMeshSubset&);

        //- Disallow default bitwise assignment
        void operator=(const newFvMeshSubset&);


        //- Check if subset has been performed
        bool checkCellSubset() const;

        //- Mark points in Map
        static void markPoints(const labelList&, Map<label>&);

        //- Mark points (with 0) in labelList
        static void markPoints(const labelList&, labelList&);

        //- Adapt nCellsUsingFace for coupled faces becoming 'uncoupled'.
        void doCoupledPatches
        (
            const bool syncPar,
            labelList& nCellsUsingFace
        ) const;

        //- Subset of subset
        static labelList subset
        (
            const label nElems,
            const labelList& selectedElements,
            const labelList& subsetMap
        );

        //- Create zones for submesh
        void subsetZones();

        //- Make dictionaries in case they are not present
        //  Creates fvSchemes and fvSolution copying from base
        void makeFvDictionaries();


public:

    // Constructors

        //- Construct given a mesh
        newFvMeshSubset(const IOobject&, const fvMesh&);


    // Destructor

        ~newFvMeshSubset();


    // Member Functions

        // Edit

            //- Set the subset. Create "oldInternalFaces" patch for exposed
            //  internal faces (patchID==-1) or use supplied patch.
            //  Does not handle coupled patches correctly if only one side
            //  gets deleted.
            void setCellSubset
            (
                const labelHashSet& globalCellMap,
                const label patchID = -1
            );

            //- Set the subset from all cells with region == currentRegion.
            //  Create "oldInternalFaces" patch for exposed
            //  internal faces (patchID==-1) or use supplied patch.
            //  Handles coupled patches by if nessecary making coupled patch
            //  face part of patchID (so uncoupled)
            void setLargeCellSubset
            (
                const labelList& region,
                const label currentRegion,
                const label patchID = -1,
                const bool syncCouples = true
            );

            //- setLargeCellSubset but with labelHashSet.
            void setLargeCellSubset
            (
                const labelHashSet& globalCellMap,
                const label patchID = -1,
                const bool syncPar = true
            );


        // Access

            //- Return reference to base mesh
            const fvMesh& baseMesh() const
            {
                return baseMesh_;
            }

            //- Return reference to subset mesh
            const fvMesh& subMesh() const;

            fvMesh& subMesh();

            //- Return reference to demand-driven subset pointMesh
            const pointMesh& subPointMesh() const;

            pointMesh& subPointMesh();

            //- Return point map
            const labelList& pointMap() const;

            //- Return face map
            const labelList& faceMap() const;

            //- Return cell map
            const labelList& cellMap() const;

            //- Return patch map
            const labelList& patchMap() const;


        // Field mapping

            //- Map volume field given map addressing
            //  - it does not patch up the exposed patches with surface fields
            //  - it does not handle patch type creation correctly
            //  Please use only in dummy mapping operations
            //  HJ, 31/Jul/2009
            template<class Type>
            static tmp<GeometricField<Type, fvPatchField, volMesh> >
            meshToMesh
            (
                const GeometricField<Type, fvPatchField, volMesh>&,
                const fvMesh& sMesh,
                const labelList& patchMap,
                const labelList& cellMap,
                const labelList& faceMap
            );

            //- Map volume field
            template<class Type>
            tmp<GeometricField<Type, fvPatchField, volMesh> >
            interpolate
            (
                const GeometricField<Type, fvPatchField, volMesh>&
            ) const;

            //- Map surface field
            template<class Type>
            static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
            meshToMesh
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&,
                const fvMesh& sMesh,
                const labelList& patchMap,
                const labelList& faceMap
            );

            template<class Type>
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
            interpolate
            (
                const GeometricField<Type, fvsPatchField, surfaceMesh>&
            ) const;

            //- Map point field
            template<class Type>
            static tmp<GeometricField<Type, pointPatchField, pointMesh> >
            meshToMesh
            (
                const GeometricField<Type, pointPatchField, pointMesh>&,
                const pointMesh& sMesh,
                const labelList& patchMap,
                const labelList& pointMap
            );

            template<class Type>
            tmp<GeometricField<Type, pointPatchField, pointMesh> >
            interpolate
            (
                const GeometricField<Type, pointPatchField, pointMesh>&
            ) const;


        //- Write data
        virtual bool writeData(Foam::Ostream&) const
        {
            return true;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "newFvMeshSubsetInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
