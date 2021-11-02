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
#include "immersedBoundaryFvPatchFieldsFwd.H"
#include "extendedNLevelGlobalCellToCellStencils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class immersedMeshMapper;

/*---------------------------------------------------------------------------*\
                       Class immersedBoundaryObject Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryObject
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

        //- Name of stl file
        const polyMesh& pMesh_;

        //- Reference to cell neighbours
        const extendedNLevelCPCCellToCellStencil& cellNeighbours_;

        //- Dictionary
        const dictionary& dict_;

        //- Used geometric dimensions
        vector geometricD_;

        //- Immersed shape
        autoPtr<immersedShape> shape_;

        //- What to sample
        const interpolationMethod mode_;

        //- Shape of immersed object
        standAlonePatch& patch_;

        //- Cells on the immersed boundary
        mutable labelList* boundaryCellsPtr_;

        //- Cells within the immersed boundary
        mutable labelList* internalCellsPtr_;

        //- Cells within the immersed boundary shell
        mutable labelList* shellCellsPtr_;

        //- Closest cells to the faces on the inside of the object
        mutable labelList patchInternalCells_;

        //- Closest cells to the faces on the outside of the object
        mutable labelList patchExternalCells_;

        //- Internal cell centres
        mutable Field<vector> internalC_;

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

        //- Nearest neighbours of immersed mesh
        labelListList nearestNeighbours_;

        //- Mass of object
        scalar mass_;

        //- Forces acting on the object
        vectorField force_;

        //- Effective force
        vector forceEff_;

        //- Effective moment
        vector momentEff_;

        //- Force due to external motion (i.e. fluid and gravity)
        vector forceExt_;

        //- Moment due to external motion (i.e. fluid and gravity)
        vector momentExt_;

        //- Gravitational acceleration
        vector g_;

        //- Switch to turn reporting of motion data on and off
        Switch report_;

        // Boundary conditions

            //- Scalar boundary conditions
            PtrListDictionary<immersedBoundaryScalarPatchField> scalarBoundaries_;

            //- Vector boundary conditions
            PtrListDictionary<immersedBoundaryVectorPatchField> vectorBoundaries_;



    // Protected Member Functions

        //- Calculate mapping
        void calcMapping() const;

        //- Update weights
        void updateWeights() const;


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
                const polyMesh& pMesh,
                const dictionary& dict,
                const dictionary& stateDict
            ),
            (pMesh, dict, stateDict)
        );

    // Constructors

        //- Construct from dictionary
        immersedBoundaryObject
        (
            const polyMesh&,
            const dictionary&,
            const dictionary&
        );


    // Selectors

        static autoPtr<immersedBoundaryObject> New
        (
            const polyMesh& pMesh,
            const dictionary& dict,
            const dictionary& stateDict
        );


    //- Destructor
    virtual ~immersedBoundaryObject();


    // Member Functions

        //- Initialize
        virtual void initialize();

        //- Clear mapping
        void clearOut();

        //- Return the global polyMesh
        const polyMesh& pMesh() const;

        //- Is the object moving
        virtual bool moving() const = 0;

        //- Does the object have a temperature dependency
        virtual bool temperatureDependent() const = 0;

        //- Return a pointer to the immersedMesh if it exists
        virtual const fvMesh& immersedMesh() const;

        //- Return a pointer to the immersedMeshMapper if it exists
        virtual const immersedMeshMapper* mapper() const;

        // Update state

            //- Update objects
            virtual bool update
            (
                bool firstIter,
                const vector& fGlobal,
                const vector& tauGlobal,
                scalar deltaT,
                scalar deltaT0
            )
            {
                return false;
            }

            //- Update
            virtual void update()
            {}

            //- Solve
            virtual void solve()
            {}

            //- Report the status of the motion
            virtual void status() const;

            //- Apply the restraints to the object
            virtual void applyRestraints() = 0;

            //- Update and relax accelerations from the force and torque
            virtual void updateAcceleration
            (
                const vector& fGlobal,
                const vector& tauGlobal
            ) = 0;

            virtual void newTime() = 0;


        // Transformations

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

            //- Apply rotation tensors to Q0 for the given torque (pi) and deltaT
            //  and return the rotated Q and pi as a tuple
            virtual Tuple2<tensor, vector> rotate
            (
                const tensor& Q0,
                const vector& pi,
                const scalar deltaT
            ) const = 0;


        // Access IBM mapping variable

            //- Return the name of the object
            inline word keyword() const;

            //- Return the name of the object
            inline const word& name() const;

            //- Return the report Switch
            inline bool report() const;

            //- Return the number of faces
            inline label nFaces() const;

            //- Return the triSurfaceMesh
            inline const standAlonePatch& patch() const;

            //- Return face centres
            inline const vectorField& faceCentres() const;

            //- Return points
            inline const pointField& points() const;

            //- Return old points
            inline const pointField& faceCentresOld() const;

            //- Return face areas * face normals
            inline const vectorField& Sf() const;

            //- Return face areas
            inline tmp<scalarField> magSf() const;

            //- Return face areas
            inline const scalarField& deltaCoeffs() const;

            //- Return the internal cells
            inline const labelList& internalCells() const;

            //- Return the face-to-cell mapping
            inline const labelList& patchInternalCells() const;

            //- Return the face-to-cell mapping
            inline const labelList& patchExternalCells() const;

            //- Return the internal cell centres
            inline const vectorField& internalC() const;

            //- Return the boundary cells
            inline const labelList& boundaryCells() const;

            //- Return the internal cells
            inline const labelList& shellCells() const;

            //- Return cells used to interpolate
            inline const labelListList& interpToCells() const;

            //- Return weights used to interpolate
            inline const scalarListList& interpToWeights() const;

            //- Return the shape
            inline const immersedShape& shape() const;

            //- Return the patch name if an immersedFvMesh is used
            virtual word patchName() const
            {
                NotImplemented;
                return word::null;
            }


        // Access physical properties

            //- Return the mass of the object
            inline scalar mass() const;

            //- Return the current centre of rotation
            virtual point centreOfRotation() const = 0;

            //- Return the current centre of mass
            inline point initialCentreOfMass() const;

            //- Return the current centre of mass
            virtual point centreOfMass() const = 0;

            //- Return the translational constraint tensor
            virtual tensor tConstraints() const = 0;

            //- Return the rotational constraint tensor
            virtual tensor rConstraints() const = 0;

            //- Return the orientation tensor, Q.
            //  globalVector = Q & bodyLocalVector
            //  bodyLocalVector = Q.T() & globalVector
            virtual tensor orientation() const = 0;

            //- Return the velocity of the center of rotation
            virtual vector v() const = 0;

            //- Return the rotated inertia tensor
            inline tensor invMomentOfInertia() const;

            //- Return the rotated inertia tensor
            virtual diagTensor momentOfInertia() const = 0;

            //- Return the velocity at a list of position
            virtual vector velocity(const point& pt) const = 0;

            //- Return the velocity at a list of position
            virtual tmp<vectorField> velocity(const pointField& pt) const = 0;

            //- Return the velocity at the faces using the previous
            //  face centres
            virtual tmp<vectorField> velocity() const;

            //- Return the velocity at the faces using the previous
            //  face centres
            virtual vector velocity(const label) const;

            //- Return the angular velocity in the global frame
            virtual vector omega() const = 0;

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

            //- Return the gravitation acceleration vector
            inline const vector& g() const;

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
            void setInternal(Field<Type>& vf, const Type& val) const;

            //- Set internal cells to the given value
            template<class Type, class CombineOp>
            void setInternal
            (
                Field<Type>& vf,
                const Type& val,
                const CombineOp& cop
            ) const;

            //- Set internal cells to the given value
            template<class Type>
            void setInternal(Field<Type>& vf, const Field<Type>& vals) const;

            //- Set internal cells to the given value
            template<class Type, class CombineOp>
            void setInternal
            (
                Field<Type>& vf,
                const Field<Type>& vals,
                const CombineOp& cop
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
                const word& name,
                GeometricField<Type, fvPatchField, volMesh>& F,
                const volScalarField& alphaRho,
                const GeometricField<Type, fvPatchField, volMesh>& alphaRhoFOld,
                const GeometricField<Type, fvPatchField, volMesh>& RHS,
                const dimensionedScalar& dT
            ) const;

            //- Add boundaries
            template<class Type>
            void addField
            (
                GeometricField<Type, fvPatchField, volMesh>&,
                const dictionary&
            );


        // Edit

            //- Move points
            virtual void movePoints()
            {}


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
