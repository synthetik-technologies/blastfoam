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
    Foam::newLeastSquaresVolPointInterpolation

Description
    Foam::newLeastSquaresVolPointInterpolation

SourceFiles
    newLeastSquaresVolPointInterpolation.C
    newLeastSquaresVolPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef newLeastSquaresVolPointInterpolation_H
#define newLeastSquaresVolPointInterpolation_H

#include "MeshObject.H"
#include "PtrList.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "scalarList.H"
#include "className.H"
#include "scalarMatrices.H"
#include "FieldField.H"
#include "Map.H"
#include "polyPatch.H"
#include "Tuple2.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
           Class newLeastSquaresVolPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

typedef labelUList unallocLabelList;

class newLeastSquaresVolPointInterpolation
:
    public MeshObject
    <
        fvMesh, UpdateableMeshObject, newLeastSquaresVolPointInterpolation
    >
{
    // Private data

        //- Point boundary face addressig
        mutable autoPtr<labelListList> pointBndFacesPtr_;

        //- Point cyclic boundary face addressig
        mutable autoPtr<labelListList> pointCyclicFacesPtr_;

        //- Point boundary face addressig
        mutable autoPtr<labelListList> pointCyclicBndFacesPtr_;

        //- Point cyclic ggi boundary face addressig
        mutable autoPtr<labelListList> pointCyclicGgiFacesPtr_;

        //- Point cyclic ggi boundary face addressig
        mutable autoPtr<labelListList> pointCyclicGgiBndFacesPtr_;

        //- Point boundary face addressig
        mutable autoPtr<labelListList> pointProcFacesPtr_;

        //- Axis edges for axisymmetric cases
        mutable autoPtr<labelList> axisEdgesPtr_;

        //- Point ngb processor boundary face centres
        mutable autoPtr<Map<labelList> > pointAxisEdgesPtr_;

        //- Global point ngb processor boundary face centres
        mutable autoPtr<Map<Field<vector> > > globalPointNgbProcBndFaceCentresPtr_;

        //- Global point ngb processor cell centres
        mutable autoPtr<Map<Field<vector> > > globalPointNgbProcCellCentresPtr_;

        //- Locel boundary faces needed by neighbour processors
        mutable autoPtr<labelListList> procBndFacesPtr_;

        //- Centres of boundary faces from neigbour processores
        mutable autoPtr<FieldField<Field, vector> > procBndFaceCentresPtr_;

        //- Point ngb processor boundary faces addressing
        mutable autoPtr<List<List<labelPair> > > pointProcBndFacesPtr_;

        //- Locel processor cells needed by neighbour processors
        mutable autoPtr<labelListList> procCellsPtr_;

        //- Point ngb processor cell addressing
        mutable autoPtr<Map<List<labelPair> > > pointProcCellsPtr_;

        //- Centres of cells from neigbour processores
        mutable autoPtr<FieldField<Field, vector> > procCellCentresPtr_;

        //- Weighting factors
        mutable autoPtr<FieldField<Field, scalar> > weightsPtr_;

        //- Local origin
        mutable autoPtr<Field<vector> > originsPtr_;

        //- Mirror plane normals and transformation tensors
        mutable autoPtr<List<Tuple2<vector, tensor> > > mirrorPlaneTransformationPtr_;

        //- Interpolation matrices
        mutable PtrList<scalarRectangularMatrix> invLsMatrices_;

        //- Interpolatin stencil size
        mutable autoPtr<scalarField> refLPtr_;

        //- Switch to check that there are processor boundaries. To run in
        //  parallel, a decomposition like 'none' needs to be used where all
        // the cells are on one processor so there are no procesosr boundaries
        Switch processorBoundariesExist_;

    // Private member functions

        //- Make point boundary face addressing
        void makePointFaces() const;

        //- Make axis edges
        void makeAxisEdges() const;

        //- Make point axis edges addressing
        void makePointAxisEdges() const;

//         //- Make point ngb processor boundary face centres
//         void makePointNgbProcBndFaceCentres() const;

        //- Make global point ngb processor cell centres
        void makeGlobalPointNgbProcBndFaceCentres() const;

        //- Make global point ngb processor cell centres
        void makeGlobalPointNgbProcCellCentres() const;

        //- Make local boundary faces needed by neighbour processors
        void makeProcBndFaces() const;

        //- Make boundary face centres from neighbour processors
        void makeProcBndFaceCentres() const;

        //- Make cells needed by neighbour processors
        void makeProcCells() const;

        //- Make cell centres from neighbour processors
        void makeProcCellCentres() const;

        //- Make weights
        void makeWeights() const;

        //- Make local origins
        void makeOrigins() const;

        //- Make interpolation matrix
        void makeInvLsMatrices() const;

        //- Make mirror plane noramals and transformation tensors
        void makeMirrorPlaneTransformation() const;

        //- foam-extend-4.0 implementation of Householder inverse
        tensor hinv(const tensor& t) const;

        //- Disallow default bitwise copy construct
        newLeastSquaresVolPointInterpolation
        (
            const newLeastSquaresVolPointInterpolation&
        );

        //- Disallow default bitwise assignment
        void operator=(const newLeastSquaresVolPointInterpolation&);

public:

    // Declare name of the class and its debug switch
    TypeName("newLeastSquaresVolPointInterpolation");


    // Constructors

        //- Constructor given fvMesh.  pointMesh will be created or
        //  looked up from objectRegistry
        explicit newLeastSquaresVolPointInterpolation(const fvMesh&);


    // Destructor

        virtual ~newLeastSquaresVolPointInterpolation();


    // Member functions

        // Access

            //- Get point boundary face addressing
            const labelListList& pointBndFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicBndFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicGgiFaces() const;

            //- Get point boundary face addressing
            const labelListList& pointCyclicGgiBndFaces() const;

            //- Get point processor boundary face addressing
            const labelListList& pointProcFaces() const;

            //- Get axis edges
            const labelList& axisEdges() const;

            //- Get point axis edges addressing
            const Map<labelList>& pointAxisEdges() const;

//             //- Get point processor boundary face addressing
//             const Map<Field<vector> >& pointNgbProcBndFaceCentres() const;

            //- Get global point ngb processor boundary face centres
            const Map<Field<vector> >&
            globalPointNgbProcBndFaceCentres() const;

            //- Get global point ngb processor cell centres
            const Map<Field<vector> >& globalPointNgbProcCellCentres() const;

            //- Get list of boundary faces needed by neighbour processores
            const labelListList& procBndFaces() const;

            //- Get centres of boundary faces from neighbour processors
            const FieldField<Field, vector>& procBndFaceCentres() const;

            //- Get neighbor processor cell point addresing
            const List<List<labelPair> >& pointProcBndFaces() const;
//             const Map<List<labelPair> >& pointProcBndFaces() const;

            //- Get list of cells needed by neighbour processores
            const labelListList& procCells() const;

            //- Get neighbor processor cell point addresing
            const Map<List<labelPair> >& pointProcCells() const;

            //- Get centres of cells from neighbour processors
            const FieldField<Field, vector>& procCellCentres() const;

            //- Get local origins
            const FieldField<Field, scalar>& weights() const;

            //- Get local origins
            const Field<vector>& origins() const;

            //- Get interpolation matrices
            const PtrList<scalarRectangularMatrix>& invLsMatrices() const;

            //- Get interpolation stencil size
            const scalarField& refL() const;

            //- Get mirror plane normals
            const List<Tuple2<vector, tensor> >&
            mirrorPlaneTransformation() const;
//             const Map<Tuple2<vector, tensor> >&
//             mirrorPlaneTransformation() const;

        // Edit

            //- Correct weighting factors for moving mesh.
            //  Updated for MeshObject.  HJ, 30/Aug/2010
            virtual bool movePoints();

            //- Update mesh topology using the morph engine
            //  Updated for MeshObject.  HJ, 30/Aug/2010
            virtual void updateMesh(const mapPolyMesh&);


    // Interpolation functions

        //- Interpolate from volField to pointField
        //  using least squares linear fit
        template<class Type>
        void interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&,
            GeometricField<Type, pointPatchField, pointMesh>&
        ) const;

        //- Interpolate from volField to pointField
        //  using least squares linear fit
        template<class Type>
        tmp<Field<Type> > interpolate
        (
            const polyPatch& patch,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Interpolate from volField to pointField
        //  using least squares linear fit
        template<class Type>
        Type interpolate
        (
            const label pointID,
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;

        //- Get neighbour processor boundary face data
        template<class Type>
        void pointNgbProcBndFaceFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Map<Field<Type> >& fieldData
        ) const;

        //- Get neighbour processor cell data
        template<class Type>
        void globalPointNgbProcBndFaceFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Map<Field<Type> >& fieldData
        ) const;

        //- Get neighbour processor cell data
        template<class Type>
        void globalPointNgbProcCellFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf,
            Map<Field<Type> >& fieldData
        ) const;

        //- Send and receive field data from neighbour processors
        template<class Type>
        tmp<FieldField<Field, Type> > procCellsFieldData
        (
            const Field<Type>& psi
        ) const;

        //- Send and receive field data from neighbour processors
        template<class Type>
        tmp<FieldField<Field, Type> > procBndFacesFieldData
        (
            const GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "newLeastSquaresVolPointInterpolate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
