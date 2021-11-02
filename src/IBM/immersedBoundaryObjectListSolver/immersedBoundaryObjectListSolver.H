/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2019 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::immersedBoundaryObjectListSolver

Description
    6-DoF solid-body mesh motion solver for an fvMesh.

    Applies SLERP interpolation of movement as function of distance to
    the object surface.

SourceFiles
    immersedBoundaryObjectListSolver.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryObjectListSolver_H
#define immersedBoundaryObjectListSolver_H

#include "immersedBoundaryObject.H"
#include "volFields.H"
#include "MeshObject.H"
#include "PtrListDictionary.H"
#include "collisionModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class immersedBoundaryObjectListSolver;

typedef MeshObject
<
    fvMesh,
    UpdateableMeshObject,
    immersedBoundaryObjectListSolver
> ImmersedBoundaryObjectListSolverObject;
/*---------------------------------------------------------------------------*\
        Class immersedBoundaryObjectListSolver Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryObjectListSolver
:
    public ImmersedBoundaryObjectListSolverObject
{
    friend immersedBoundaryObject;

    // Private Data

        //- Const reference to mesh
        const fvMesh& mesh_;

        //- Maximum number of times to expand cell neighbours
        label maxLevel_ = 2;

        //- Immersed boundary dictionary
        IOdictionary ibmDict_;

        //- Six DoF motion object
        PtrListDictionary<immersedBoundaryObject> objects_;


        //- Switch for test-mode in which only the
        //  gravitational body-force is applied
        Switch test_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Cell Types
        volScalarField cellTypes_;

        //- Cell Types
        autoPtr<volScalarField> objectIDPtr_;

        //- Are any objects temperature based
        bool thermalForcingNeeded_;

        //- Collision model
        autoPtr<collisionModel> collision_;

        //- Return the effective viscous stress (laminar + turbulent).
        tmp<volSymmTensorField> devTau(const word& phaseName) const;

public:

    //- Runtime type information
    TypeName("immersedBoundaryObjectListSolver");


    // Constructors

        //- Construct from polyMesh and dictionary
        immersedBoundaryObjectListSolver
        (
            const fvMesh&
        );

        //- Disallow default bitwise copy construction
        immersedBoundaryObjectListSolver
        (
            const immersedBoundaryObjectListSolver&
        ) = delete;


    //- Destructor
    virtual ~immersedBoundaryObjectListSolver();


    // Member Functions

        //- Return the six DoF motion object
        const PtrListDictionary<immersedBoundaryObject>& objects() const;

        //- Return if thermal forcing is needed
        bool thermalForcingNeeded() const
        {
            return thermalForcingNeeded_;
        }

        //- Solve for motion
        virtual void solve();

        //- Move points
        virtual bool movePoints();

        //- Update
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Add fields to correct
        template<class Type>
        void add(GeometricField<Type, fvPatchField, volMesh>& f);

        //- Update boundaries
        void correctBoundaryConditions();

        //- Return momentum forcing
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> forcing
        (
            const GeometricField<Type, fvPatchField, volMesh>& F,
            const volScalarField& alphaRho,
            const GeometricField<Type, fvPatchField, volMesh>& alphaRhoFOld,
            const GeometricField<Type, fvPatchField, volMesh>& RHS,
            const dimensionedScalar& dT
        ) const;

        //- Set an internal field to a given value
        template<class Type>
        void setInternal
        (
            Field<Type>& f,
            const Type& v
        )
        {
            forAll(objects_, i)
            {
                objects_[i].setInternal(f, v);
            }
        }

        //- Return the Courant number (max of all objects)
        scalar CoNum() const;

        //- Return the maximum Courant number
        scalar maxCoNum() const;

        //- Write motion state information for restart
        virtual bool write(const bool write = true) const;
        virtual bool writeData(Ostream&) const
        {
            return write(true);
        }
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp,
            const bool write
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const immersedBoundaryObjectListSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "immersedBoundaryObjectListSolverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
