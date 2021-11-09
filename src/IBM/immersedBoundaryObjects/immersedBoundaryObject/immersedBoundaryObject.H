/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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
    Foam::immersedBoundaryObject

Description


SourceFiles
    immersedBoundaryObject.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryObject_H
#define immersedBoundaryObject_H

#include "polyMesh.H"
#include "immersedShape.H"
#include "PtrListDictionary.H"
#include "volFields.H"
#include "primitiveFieldsFwd.H"
#include "extendedNLevelGlobalCellToCellStencils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class immersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryObject
:
    public standAlonePatch
{
    // Private member functions

        //- Return the delta function
        scalar delta
        (
            const point& p,
            const point& c,
            const scalar& h
        ) const;


public:

    // Type enumerations

        //- Mesh items to sample
        enum interpolationMethod
        {
            INVERSE_DISTANCE
        };

        static const NamedEnum<interpolationMethod, 1>
            interpolationMethodNames_;


protected:

    // Protected data

        //- Reference to the associated polyPatch
        const polyPatch& patch_;

        //- Name of stl file
        const polyMesh& pMesh_;

        //- Dictionary
        const dictionary& dict_;

        //- What to sample
        const interpolationMethod mode_;

        //- Immersed shape
        mutable autoPtr<immersedShape> shape_;

        //- Cells on the immersed boundary
        mutable labelList* boundaryCellsPtr_;

        //- Cells within the immersed boundary
        mutable labelList* allInternalCellsPtr_;

        //- Cells within the immersed boundary
        mutable labelList* internalCellsPtr_;

        //- Cells within the immersed boundary shell
        mutable labelList* shellCellsPtr_;

        //- Closest cells to the faces on the inside of the object
        mutable labelList patchInternalCells_;

        //- Closest cells to the faces on the outside of the object
        mutable labelList patchExternalCells_;

        //- Closest cells to the object
        //  patchExternal cells unless one does not exist
        mutable  labelList patchCells_;

        //- Map from patch cells to the face index
        mutable  labelListList patchMap_;

        //- Cells containing each surface point
        mutable List<labelList> interpToCells_;

        //- Interpolation weights
        mutable List<scalarList> interpToWeights_;

        //- Cells forcing is applied to (for each surface point)
        mutable List<labelList> interpFromPoints_;

        //- Forcing weights
        mutable List<scalarList> interpFromWeights_;

        //- Delta coefficients
        //  Inverse distance from closest cell to face centre
        mutable scalarField deltaCoeffs_;

        mutable scalarField* magSfPtr_;

        //- Nearest neighbours of immersed mesh
        mutable labelListList nearestNeighbours_;

        //- Forces acting on the object
        mutable vectorField force_;

        //- Effective force
        vector forceEff_;

        //- Effective moment
        vector momentEff_;

        //- Force due to external motion (i.e. fluid and gravity)
        vector forceExt_;

        //- Moment due to external motion (i.e. fluid and gravity)
        vector momentExt_;



    // Protected Member Functions

        //- Initialize
        virtual void initialize() const;

        //- Calculate mapping
        void calcMapping() const;

        //- Update weights
//         void updateWeights() const;


public:

    //- Runtime type information
    ClassName("immersedBoundaryObject");

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            immersedBoundaryObject,
            dictionary,
            (
                const polyPatch& patch,
                const dictionary& dict,
                const dictionary& stateDict
            ),
            (patch, dict, stateDict)
        );

    // Constructors

        //- Construct from dictionary
        immersedBoundaryObject
        (
            const polyPatch&,
            const dictionary&,
            const dictionary&
        );


    // Selectors

        static autoPtr<immersedBoundaryObject> New
        (
            const polyPatch& patch,
            const dictionary& dict,
            const dictionary& stateDict
        );


    //- Destructor
    virtual ~immersedBoundaryObject();


    // Member Functions

        //- Clear mapping
        void clearOut();

        //- Return the global polyMesh
        inline const polyMesh& pMesh() const;

        //- Return the polyPatch
        inline const polyPatch& patch() const;

        //- Return the patch index
        inline label index() const;

        //- Is the object moving
        virtual bool moving() const = 0;

        //- Does the object have a temperature dependency
        virtual bool temperatureDependent() const = 0;

        //- Update
        virtual void update()
        {}

        //- Solve
        virtual void solve()
        {}

        //- Report the status of the motion
        virtual void status(const bool print = true) const;


        // Transformations

            //- Return the orientation tensor, Q.
            //  globalVector = Q & bodyLocalVector
            //  bodyLocalVector = Q.T() & globalVector
            virtual tensor orientation() const = 0;

            //- Transform the given initial state point by the current motion
            //  state
            virtual point transform(const point&) const = 0;

            //- Transform the given initial state pointField by the current
            //  motion state
            virtual tmp<pointField> transform(const pointField&) const = 0;

            //- Inverse of transform
            virtual point inverseTransform(const point&) const = 0;

            //- Inverse of transform
            virtual tmp<pointField> inverseTransform(const pointField&) const = 0;


        // Access IBM mapping variable

            //- Return the name of the object
            inline const word& name() const;

            //- Return the name of the object
            inline const word& keyword() const;

            //- Return face areas
            inline const scalarField& deltaCoeffs() const;

            //- Return all the internal cells
            inline const labelList& allInternalCells() const;

            //- Return the internal cells
            inline const labelList& internalCells() const;

            //- Return the face-to-cell mapping
            inline const labelList& patchInternalCells() const;

            //- Return the face-to-cell mapping
            inline const labelList& patchExternalCells() const;

            //- Return the face-to-cell mapping
            inline const labelList& patchCells() const;

            //- Return the map from the patch cells to the face
            inline const labelListList& patchMap() const;

            //- Return the boundary cells
            inline const labelList& boundaryCells() const;

            //- Return the internal cells
            inline const labelList& shellCells() const;

            //- Return cells used to interpolate
            inline const labelListList& interpToCells() const;

            //- Return weights used to interpolate
            inline const scalarListList& interpToWeights() const;

            //- Return the number of faces
            inline label nFaces() const;

            //- Return the face area vector
            inline const vectorField& Sf() const;

            //- Return the face area
            inline const scalarField& magSf() const;

            //- Return the face normals
            inline tmp<vectorField> nf() const;

            //- Return the shape
            inline const immersedShape& shape() const;

            //- Store the motion state at the beginning of the time-step
            virtual void newTime()
            {}


        // Access physical properties

            //- Return the mass of the object
            virtual scalar mass() const = 0;

            //- Return the current centre
            virtual point centre() const = 0;

            //- Return the rotated inertia tensor
            inline tensor invMomentOfInertia() const;

            //- Return the rotated inertia tensor
            virtual diagTensor momentOfInertia() const = 0;

            //- Return the velocity at a location
            virtual vector v(const point& pt) const = 0;

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity(const pointField& pt) const = 0;

            //- Return the velocity at the faces using the previous
            //  face centres
            virtual tmp<vectorField> velocity() const = 0;

            //- Return the force
            inline const vectorField& force() const;

            //- Access the force
            inline vectorField& force();

            //- Return the effective force
            inline const vector& forceEff() const;

            //- Access the effective force
            inline vector& forceEff();

            //- Return the effective moment
            inline const vector& momentEff() const;

            //- Access the effective moment
            inline vector& momentEff();

            //- Return the external force
            inline const vector& forceExt() const;

            //- Access the external force
            inline vector& forceExt();

            //- Return the external moment
            inline const vector& momentExt() const;

            //- Access the external moment
            inline vector& momentExt();

            //- Return the current Courant number
            virtual scalar CoNum() const
            {
                return 0.0;
            }

            //- Return the max Courant number
            virtual scalar maxCoNum() const
            {
                return great;
            }


        // Interpolation and forcing

            //- Interpolate using weights
            template<class Type>
            tmp<Field<Type>> interpolateTo(const Field<Type>&) const;

            //- Boundary values using the closest cell (inside)
            template<class Type>
            tmp<Field<Type>> patchInternalField(const Field<Type>&) const;

            //- Boundary values using the closest cell (outside)
            template<class Type>
            tmp<Field<Type>> patchExternalField(const Field<Type>&) const;

            //- Boundary values using the closest cell (outside)
            template<class Type>
            tmp<Field<Type>> patchField(const Field<Type>&) const;

            //- Interpolate using weights
            template<class Type>
            tmp<Field<Type>> boundaryValues(const Field<Type>&) const;

            //- Apply forcing to volume field given the surface force
            template<class Type>
            void interpolateFrom
            (
                const Field<Type>& sf,
                Field<Type>& vf,
                const bool clear = false
            ) const;

            //- Set internal cells to the given value
            template<class Type>
            void setInternal
            (
                Field<Type>& vf,
                const Type& val,
                const bool all = false
            ) const;

            //- Set internal cells to the given value
            template<class Type, class CombineOp>
            void setInternal
            (
                Field<Type>& vf,
                const Type& val,
                const CombineOp& cop,
                const bool all = false
            ) const;

            //- Set internal cells to the given value
            template<class Type>
            void setInternal
            (
                Field<Type>& vf,
                const Field<Type>& vals,
                const bool all = false
            ) const;

            //- Set internal cells to the given value
            template<class Type, class CombineOp>
            void setInternal
            (
                Field<Type>& vf,
                const Field<Type>& vals,
                const CombineOp& cop,
                const bool all = false
            ) const;

            //- Set forcing cells to the given value
            template<class Type>
            void setShell(Field<Type>& vf, const Type& val) const;

            //- Set forcing cells to the given value
            template<class Type, class CombineOp>
            void setShell
            (
                Field<Type>& vf,
                const Type& val,
                const CombineOp& cop
            ) const;

            //- Set boundary cells to the given value
            template<class Type>
            void setBoundary(Field<Type>& vf, const Type& val) const;

            //- Set boundary cells to the given value
            template<class Type, class CombineOp>
            void setBoundary
            (
                Field<Type>& vf,
                const Type& val,
                const CombineOp& cop
            ) const;

            //- Return if points are found in the object
            labelList calcInside(const pointField& points) const;

            //- Return the bounding box
            inline const boundBox& bounds() const;

            //- Set values inside the object using the boundary conditions
            void setValues();

            //- Set boundary vales
            template<class Type>
            void setValues(DimensionedField<Type, volMesh>& field);

            //- Return momentum forcing
            template<class Type>
            void addForcing
            (
                Field<Type>& F,
                const Field<Type>& pF,
                const Field<scalar>& alphaRho,
                const Field<Type>& alphaRhoFOld,
                const Field<Type>& RHS,
                const scalar dT
            ) const;


        // Edit

            using standAlonePatch::movePoints;

            //- Move points
            virtual void movePoints();


        // I/O

            //- Write
            virtual void write(Ostream&) const;

            virtual void write(dictionary&) const;

            //- Read coefficients dictionary and update system parameters,
            //  constraints and restraints but not the current state
            virtual bool read(const dictionary& dict);
};

//- Add to list if the entry is not found
bool addIfUnique
(
    labelList& l,
    const label li,
    const label entry
);

//- Add to list if the entry is not found and increment
bool addIfUniqueAndIncrement
(
    labelList& l,
    label& li,
    const label entry
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "immersedBoundaryObjectTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "immersedBoundaryObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
