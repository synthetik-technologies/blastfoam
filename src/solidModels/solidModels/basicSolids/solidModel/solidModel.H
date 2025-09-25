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
    solidModel

Description
    Virtual base class for solid mechanics models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    solidModel.C
    solidModelTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef solidModel_H
#define solidModel_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "scalarIOField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "thermalModel.H"
#include "mechanicalModel.H"
#include "nonLinearGeometry.H"
#include "uniformDimensionedFields.H"
#include "momentumStabilisation.H"
#include "setCellDisplacements.H"
#include "OFstream.H"
#include "fvMatricesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class globalPolyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class solidModel Declaration
\*---------------------------------------------------------------------------*/

class solidModel
:
//     public physicsModel,
    public IOdictionary
{
protected:
    // Protected data

        //- Mesh
        fvMesh& mesh_;

        //- Derived solidModel type
        const word type_;

        //- Thermal model
        thermalModel thermal_;

        //- Reference to globalPolyBoundaryMesh object
        globalPolyBoundaryMesh& globalPatches_;


private:
    // Private member functions

        //- Disallow default bitwise copy construct
        solidModel(const solidModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidModel&);


protected:

    // Protected member functions

        // Check wedge definition
        static void checkWedges(const volVectorField& D);

        // Return point boundary conditions
        static wordList pointDBoundaryTypes(const volVectorField& D);

        //- This function will check that the field was read from disk
        void isRequired(const volVectorField& D, const word& type) const;

        // Read from dictionary
        virtual void readDict();

        //- Return non-const reference to solidModelCoeffs dictionary
        dictionary& solidModelDict();

        //- Return non-const reference to the thermal model
        thermalModel& thermal();

        //- Return non-const reference to the density
        volScalarField& rho();


public:

    //- Runtime type information
    TypeName("solidModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidModel,
            dictionary,
            (
                fvMesh& mesh
            ),
            (mesh)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidModel,
            lagrangian,
            (
                fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        solidModel
        (
            const word& type,
            fvMesh& mesh,
            const bool isSolid = true
        );


    // Selectors

        //- Select constructed from time
        static autoPtr<solidModel> New(fvMesh& mesh);

        //- Select constructed from time (moving meshes)
        static autoPtr<solidModel> NewLU(fvMesh& mesh);


    // Destructor

        virtual ~solidModel();


    // Member Functions

        // Access

            //- Return time
            inline const Time& runTime() const
            {
                return mesh_.time();
            }

            //- Return mesh
            inline const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return non-const mesh
            inline fvMesh& mesh()
            {
                return mesh_;
            }

            //- Return point mesh
            inline const pointMesh& pMesh() const
            {
                return pointMesh::New(mesh_);
            }

            //- Return non-const point mesh
            inline pointMesh& pMesh()
            {
                return pointMesh::New(mesh_);
            }

            //- Return the global polyPatches
            inline const globalPolyBoundaryMesh& globalPatches() const
            {
                return globalPatches_;
            }

            //- Access the global polyPatches
            inline globalPolyBoundaryMesh& globalPatches()
            {
                return globalPatches_;
            }

            // Initialize global poly patch displacement
            void initialize();

            //- Return const-reference to the density
            const volScalarField& rho() const;

            //- Return const reference to solidModelCoeffs dictionary
            const dictionary& solidModelDict() const;

            //- Return const reference to the thermal model
            const thermalModel& thermal() const;

            //- Return const reference to D
            virtual const volVectorField& D() const = 0;

            //- Return const reference to U
            virtual const volVectorField& U() const = 0;

            //- Return const reference to pointD
            virtual const pointVectorField& pointD() const = 0;

            //- Return nonlinear geometry enumerator
            virtual nonLinearGeometry::nonLinearType nonLinGeom() const = 0;

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual const volVectorField& solutionD() const = 0;

             //- Each solidModel must indicate whether gradD or gradDD
            //  is the primary solution variable
            virtual const volTensorField& solutionGradD() const = 0;

            //- Return enforceLinear switch
            virtual bool enforceLinear() const = 0;

            //- Return enforceLinear switch
            virtual void enforceLinear(const bool tf) = 0;

            //- Does this model take an incremental approach?
            //  i.e. does it solve for DD instead of D?
            //  This defaults to false but can be overwritten
            virtual bool incremental() const
            {
                return false;
            }

            //- Does this model move the mesh?
            virtual bool movingMesh() const
            {
                return false;
            }

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


        // Edit

            //- Update fields
            virtual void update(const bool correctSigma = true) = 0;

            //- Check if solid model is diverging using the Jacobian field
            virtual bool checkEnforceLinear(const volScalarField& J);

            //- Check if solid model is diverging using the Jacobian field
            virtual bool checkEnforceLinear(const surfaceScalarField& J);

            //- Evolve the solid model
            virtual bool evolve() = 0;

            //- Traction boundary surface normal gradient
            //  uses tractionBoundary for boundary traction
            virtual tmp<Foam::vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Boundary traction
            virtual tmp<Foam::vectorField> nf(const fvPatch& patch) const = 0;

            //- Return the patch implicit stiffness
            virtual tmp<scalarField> impK(const fvPatch& patch) const = 0;

            //- Return the patch stress
            virtual const symmTensorField& sigma(const fvPatch& patch) const = 0;

            //- Update total accumulated fields
            virtual void updateTotalFields();

            //- Return the desired new time-step size
            virtual scalar newDeltaT() const;

            //- Move the mesh to the deformed configuration
            virtual void moveMesh
            (
                const pointField& oldPoints,
                const volVectorField& DD,
                pointVectorField& pointDD
            );

            //- Read from dictionary
            virtual bool read();

//             //- Read from dictionary
//             virtual bool readIfModified();

            //- Write fields
            virtual bool write(const bool write = true) const;

            virtual bool writeData(Ostream&) const
            {
                return true;
            }

            virtual bool writeObject
            (
                IOstream::streamFormat fmt,
                IOstream::versionNumber ver,
                IOstream::compressionType cmp,
                const bool write
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
