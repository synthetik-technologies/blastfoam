/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::mechanicalModel

Description
    The mechanical model class takes care of reading the mechanical properties
    and creating mechanical laws.

    The solidSubMeshes class takes care of multiple mechanical laws, where each
    each material region has its own mechanical law.

    Corrections are applied at bi-material interfaces to ensure continuity of
    stress without oscillations.

SourceFiles
    mechanicalModel.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef mechanicalModel_H
#define mechanicalModel_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "mechanicalLaw.H"
#include "Switch.H"
#include "volPointInterpolation.H"
#include "nonLinearGeometry.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class mechanicalModel Declaration
\*---------------------------------------------------------------------------*/

class mechanicalModel
:
    public IOdictionary,
    public PtrList<mechanicalLaw>
{
    // Private data

        //- Const reference to the mesh
        const fvMesh& mesh_;

        //- Is the current case plane stress or plain strain/3-D
        //  Todo: this should be moved to actual mechanical law
        const Switch planeStress_;

        //- Does the solidModel use an incremental approach?
        //  i.e. does the solidModel solve for DD instead of D
        const bool incremental_;

        //- Name of the cellZone for each material
        mutable wordList cellZoneNames_;

        //- Volume field to point field interpolator for the base mesh
        mutable volPointInterpolation* volToPointPtr_;

       //- The implicit stiffness surface field for Rhie-Chow correction
        mutable surfaceScalarField* impKfcorrPtr_;


    // Private Member Functions

        //- Make material sub-meshes
        void makeSolSubMeshes() const;

        //- Make base mesh vol-to-point interpolator
        void makeVolToPoint() const;

        //- Calculate the impKf field for Rhie-Chow correction
        //  This field is zero on bi-material interfaces
        void calcImpKfcorr() const;

        //- Return the implicit stiffness surface field for Rhie-Chow correction
        const surfaceScalarField& impKfcorr() const;

        //- Clear out demand driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        mechanicalModel(const mechanicalModel&);

        //- Disallow default bitwise assignment
        void operator=(const mechanicalModel&);


public:

    // Constructors

        //- Construct from an fvMesh
        mechanicalModel
        (
            const fvMesh& mesh,
            const nonLinearGeometry::nonLinearType& nonLinGeom,
            const bool incremental = false
        );


    //- Destructor
    ~mechanicalModel();


    // Member Functions

        // Access

            //- Return const access to the mesh
            const fvMesh& mesh() const;

            //- Return a const reference to the vol-to-point base mesh
            //  interpolator
            const volPointInterpolation& volToPoint() const;

            //- Return the implicit stiffness
            //  This is the diffusivity for the Laplacian term
            tmp<volScalarField> impK() const;

            //- Return the implicit stiffness surface field
            //  This is the diffusivity for the Laplacian term
            tmp<surfaceScalarField> impKf() const;

            //- Return bulk modulus
            tmp<volScalarField> bulkModulus() const;

            //- Return material residual i.e. a measured of how convergence of
            //  the material model
            scalar residual();

        // Edit

            //- Calculate the stress volField
            void correct(volSymmTensorField& sigma);

            //- Calculate the stress surfaceField
            void correct(surfaceSymmTensorField& sigma);

            //- Update the gradient of displacement vol field
            void grad
            (
                const volVectorField& D,
                volTensorField& gradD
            );

            //- Update the gradient of displacement vol field using pointD
            void grad
            (
                const volVectorField& D,
                const pointVectorField& pointD,
                volTensorField& gradD
            );

            //- Update the gradient of displacement surface field using pointD
            void grad
            (
                const volVectorField& D,
                const pointVectorField& pointD,
                surfaceTensorField& gradDf
            );

            //- Update the gradient of displacement vol and surface field using
            //  pointD
            void grad
            (
                const volVectorField& D,
                const pointVectorField& pointD,
                volTensorField& gradD,
                surfaceTensorField& gradDf
            );

            //- Interpolate D to pointD
            void interpolate
            (
                const volVectorField& D,
                pointVectorField& pointD,
                const bool useVolFieldSigma = true
            );

            //- Return explicit Rhie-Chow correction field to quell oscillations
            tmp<volVectorField> RhieChowCorrection
            (
                const volVectorField& D,
                const volTensorField& gradD,
                const surfaceScalarField& gamma
            ) const;

            //- Return explicit Rhie-Chow correction field to quell oscillations
            tmp<volVectorField> RhieChowCorrection
            (
                const volVectorField& D,
                const volTensorField& gradD
            ) const;

            //- Update total accumulated fields
            void updateTotalFields();

            //- Return the desired new time-step size
            scalar newDeltaT();

            //- Move subMeshes
            void moveSubMeshes();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
