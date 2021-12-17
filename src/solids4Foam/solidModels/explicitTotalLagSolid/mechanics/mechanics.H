/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::mechanics

Description
    Evaluate quantities widely used in continuum mechanics.

SourceFiles
    mechanics.C

\*---------------------------------------------------------------------------*/

#ifndef mechanics_H
#define mechanics_H

#include "fvc.H"
#include "operations.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mechanics
:
    public MeshObject<fvMesh, MoveableMeshObject, mechanics>
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Operations class
        const operations& ops_;

        //- Deformation tensor
        const volTensorField& F_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Relative deformation tensor
        volTensorField relF_;

        //- Cofactor of deformation
        volTensorField H_;

        //- Material normal
        surfaceVectorField N_;

        //- Spatial normal
        surfaceVectorField n_;

        //- Stabilisation matrix for linear momentum
        surfaceTensorField stabRhoU_;

        //- Stabilisation matrix for traction
        surfaceTensorField stabTraction_;

        //- Stretch
        volScalarField stretch_;

private:

    //- Disallow default bitwise copy construct
        mechanics(const mechanics&);

    //- Disallow default bitwise assignment
        void operator=(const mechanics&);

public:

    //- Runtime type information
    TypeName("mechanics");

    // Constructors

        //- Construct from mesh and dictionary
        mechanics
        (
            const volTensorField&,
            const operations&
        );


    //- Destructor
    virtual ~mechanics();


    // Member functions

        // Access

            //- Return the deformation tensor
            const volTensorField& F() const
            {
                return F_;
            }

            //- Return the Jacobian of deformation
            const volScalarField& J() const
            {
                return J_;
            }

            //- Return the cofactor of deformation
            const volTensorField& H() const
            {
                return H_;
            }

            //- Return the relative deformation tensor
            const volTensorField& relF() const
            {
                return relF_;
            }

            //- Return the mesh surface normals
            const surfaceVectorField& N() const
            {
                return N_;
            }

            //- Return the true surface normals
            const surfaceVectorField& n() const
            {
                return n_;
            }

            //- Momentum stabilization
            const surfaceTensorField& stabRhoU() const
            {
                return stabRhoU_;
            }

            //- Traction stabilization
            const surfaceTensorField& stabTraction() const
            {
                return stabTraction_;
            }

            //- Return the stretch
            const volScalarField stretch() const
            {
                return stretch_;
            }

        // Edit

            //- Computes spatial surface normals
            void correctN();

            //- Compute stabilisation
            void correct
            (
                const volScalarField&,
                const volScalarField&
            );

            virtual bool movePoints();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
