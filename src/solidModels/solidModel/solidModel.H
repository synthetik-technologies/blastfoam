/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
21-06-2021 Jeff Heylmun     | Created class for solid model
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    Foam::solidModel

Description
    Base class for solving solid displacement

SourceFiles
    solidModel.C
    newSolidModel.C

\*---------------------------------------------------------------------------*/

#ifndef solidModel_H
#define solidModel_H

#include "integrationSystem.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "solidBlastThermo.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidModel Declaration
\*---------------------------------------------------------------------------*/

class solidModel
:
    public IOdictionary
{
protected:
// Protected data

        //- Mesh
        fvMesh& mesh_;

        //- Derived solidModel type
        const word type_;

        //- Thermo dictionary
        IOdictionary thermophysicalProperties_;

        //- Thermodynamic model
        autoPtr<solidBlastThermo> thermoPtr_;

        //- Total displacement field
        volVectorField D_;

        //- Point total displacement field
        pointVectorField pointD_;

        //- Stress field
        //  This is the engineering stress for linear geometry approaches, and
        //  true (Cauchy) stress for nonlinear geometry (large strain)
        //  approaches
        volSymmTensorField sigma_;

        //- Solution standard tolerance
        const scalar solutionTol_;

        //- Maximum number of momentum correctors
        const int nCorr_;

        //- fvModels
        fvModels& fvModels_;

        //- fvConstraints
        fvConstraints& fvConstraints_;


        //- This function will check that the D field was read from disk
        virtual void DisRequired();

        //- This function will check that the DD field was read from disk
        virtual void DDisRequired();

        //- Disallow default bitwise copy construct
        solidModel(const solidModel&);

        //- Disallow default bitwise assignment
        void operator=(const solidModel&);


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
        solidModel(const word& type, fvMesh& mesh);


    // Selectors

        //- Select constructed from time
        static autoPtr<solidModel> New(fvMesh& mesh);

        //- Select constructed from time (moving meshes)
        static autoPtr<solidModel> NewLU(fvMesh& mesh);


    // Destructor

        virtual ~solidModel();


    // Member Functions

        // Access

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return non-const mesh
            fvMesh& mesh()
            {
                return mesh_;
            }

            //- Return const-reference to the density
            const solidBlastThermo& thermo() const;

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

            //- Evolve the solid model
            virtual bool evolve() = 0;

            //- Write function must be defined for regIOobjects
            virtual bool writeData(Ostream&) const
            {
                return false;
            }

        // Static functions

            //- Return a uniform or read field
            static tmp<volScalarField> uniformOrRead
            (
                const fvMesh& mesh,
                const dictionary& dict,
                dimensionedScalar val
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
