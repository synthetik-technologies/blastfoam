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
    DSolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef DSolidModel_H
#define DSolidModel_H

#include "solidModel.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "thermalModel.H"
#include "mechanicalModel.H"
#include "nonLinearGeometry.H"
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

class DSolidModel
:
    public solidModel
{
protected:
    // Protected data

        //- Mechanical model
        mechanicalModel mechanical_;

        //- Total displacement field
        volVectorField D_;

        //- Increment of displacement field
        //  DD = D - D.oldTime()
        volVectorField DD_;

        //- Velocity field
        volVectorField U_;

        //- Point total displacement field
        pointVectorField pointD_;

        //- Point increment of displacement field
        //  pointDD = pointD - pointD.oldTime()
        pointVectorField pointDD_;

        //- Gradient of total displacement
        volTensorField gradD_;

        //- Gradient of the displacement increment
        volTensorField gradDD_;

        //- Stress field
        //  This is the engineering stress for linear geometry approaches, and
        //  true (Cauchy) stress for nonlinear geometry (large strain)
        //  approaches
        volSymmTensorField sigma_;

        //- Gravitational acceleration
        uniformDimensionedVectorField g_;

        //- Stabilisation term for momentum equation
        autoPtr<momentumStabilisation> stabilisationPtr_;

        //- Optional: set displacement at internal cells
        mutable autoPtr<setCellDisplacements> setCellDispsPtr_;

        //- Enforce linear constraint
        bool enforceLinear_;


    // Protected Member Functions

        //- For calculating the position of the patch/face zone in the deformed
        //  configuration, we need to move the mesh by the point displacement
        //  field. For moving mesh approach (updated Lagrangian) we need the
        //  pointDD field, whereas for linear geometry (small strain) and total
        //  Lagrangian approaches, we need the pointD field. This function will
        //  check which field to return
        const pointVectorField& pointDorPointDD() const;

        //- Make setCellDisplacements
        void makeSetCellDisps() const;


        //- Return non-const reference to the mechanical model
        mechanicalModel& mechanical();

        //- Optional: cells can be forced to a specified displacement in the
        //  linear system
        void setCellDisps(fvVectorMatrix& DEqn);

        //- Return const reference to stabilisation
        const momentumStabilisation& stabilisation() const
        {
            return stabilisationPtr_();
        }

        //- Set the displacement fields if a velocity was read
        void displacementFromVelocity
        (
            volVectorField& disp,
            volVectorField& ddisp
        );

        //- Construct null
        DSolidModel(const word& type, fvMesh& mesh);

public:

    //- Runtime type information
    TypeName("DSolidModel");

    // Constructors

        //- Construct from components
        DSolidModel
        (
            const word& type,
            fvMesh& mesh,
            const nonLinearGeometry::nonLinearType nonlinear,
            const bool incremental,
            const bool isSolid = true
        );


    // Destructor

        virtual ~DSolidModel();


    // Member Functions

        // Access

            // Initialize global poly patch displacement
            virtual void initialize();

            //- Return non-const reference to D
            virtual volVectorField& D()
            {
                return D_;
            }

            //- Return const reference to D
            virtual const volVectorField& D() const
            {
                return D_;
            }

            //- Return non-const reference to DD
            virtual volVectorField& DD()
            {
                return DD_;
            }

            //- Return const reference to DD
            virtual const volVectorField& DD() const
            {
                return DD_;
            }

            //- Return non-const reference to U
            virtual volVectorField& U()
            {
                return U_;
            }

            //- Return const reference to U
            virtual const volVectorField& U() const
            {
                return U_;
            }

            //- Return non-const reference to pointD
            virtual pointVectorField& pointD()
            {
                return pointD_;
            }

            //- Return const reference to pointD
            virtual const pointVectorField& pointD() const
            {
                return pointD_;
            }

            //- Return non-const reference to pointDD
            virtual pointVectorField& pointDD()
            {
                return pointDD_;
            }

            //- Return const reference to pointDD
            virtual const pointVectorField& pointDD() const
            {
                return pointDD_;
            }

            //- Return non-const reference to gradD
            virtual volTensorField& gradD()
            {
                return gradD_;
            }

            //- Return const reference to gradD
            virtual const volTensorField& gradD() const
            {
                return gradD_;
            }

            //- Return non-const reference to gradDD
            virtual volTensorField& gradDD()
            {
                return gradDD_;
            }

            //- Return const reference to gradDD
            virtual const volTensorField& gradDD() const
            {
                return gradDD_;
            }

            //- Return non-const reference to sigma
            virtual volSymmTensorField& sigma()
            {
                return sigma_;
            }

            //- Return const reference to sigma
            virtual const volSymmTensorField& sigma() const
            {
                return sigma_;

            }

            //- Return the patch stress
            virtual const symmTensorField& sigma(const fvPatch& patch) const
            {
                return sigma_.boundaryField()[patch.index()];
            }

            //- Return const reference to the mechanical model
            const mechanicalModel& mechanical() const;

            //- Return const reference to the gravitational acceleration
            inline const uniformDimensionedVectorField& g() const
            {
                return g_;
            }

            //- Return enforceLinear switch
            virtual bool enforceLinear() const
            {
                return enforceLinear_;
            }

            //- Return enforceLinear switch
            virtual void enforceLinear(const bool tf)
            {
                enforceLinear_ = tf;
            }

            //- Return a reference to setCellDisplacements
            const setCellDisplacements& setCellDisps() const;


        // Edit

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
