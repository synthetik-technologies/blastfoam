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
    Foam::gradientSchemes

Description
    Define gradient computation schemes.

SourceFiles
    gradientSchemes.C

\*---------------------------------------------------------------------------*/

#ifndef gradientSchemes_H
#define gradientSchemes_H

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "operations.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class gradientSchemes
:
    public MeshObject<fvMesh, MoveableMeshObject, gradientSchemes>
{
    // Private data

        //- Mesh
        const fvMesh& mesh_;

        //- Operations
        const operations ops_;

        //- Cell owners
        const labelUList& own_;

        //- Cell neighbours
        const labelUList& nei_;

        // Cell centre coordinates
        const volVectorField& C_;

        // Face centre coordinates
        const surfaceVectorField& Cf_;

        //- Material nodal coordinates
        const pointField& points_;

        //- Inverse of distance matrix
        tensorField Ainv_;

        //- Inverse of local distance matrix
        tensorField AinvLocal_;

private:

    // Private member functions

        //- Disallow default bitwise copy construct
        gradientSchemes(const gradientSchemes&);

        //- Disallow default bitwise assignment
        void operator=(const gradientSchemes&);


public:

    //- Runtime type information
    TypeName("gradientSchemes");

    // Constructors

        //- Construct from mesh
        gradientSchemes(const fvMesh&);


    // Destructor
    virtual ~gradientSchemes();


    // Member functions

        // Edit

            //- Stabilised tensor inverse
            static tensor stabInv(const tensor&);

            //- Inverse distance matrix for gradient
            tmp<tensorField> distanceMatrix() const;

            //- Inverse distance matrix for local gradient
            tmp<tensorField> distanceMatrixLocal() const;

            //- Least square gradient for a volScalarField
            tmp<volVectorField> gradient
            (
                const GeometricField<scalar, fvPatchField, volMesh>& U
            ) const;

            //- Least square gradient for a volVectorField
            tmp<volTensorField> gradient
            (
                const volVectorField& U
            ) const;

            //- Least square gradient for a volTensorField
            void gradient
            (
                const volTensorField&,
                volTensorField&,
                volTensorField&,
                volTensorField&
            ) const;

            //- Local least square gradient for a volVectorField
            tmp<volTensorField> localGradient
            (
                const volVectorField&,
                const GeometricField<vector, fvsPatchField, surfaceMesh>&,
                const GeometricField<vector, pointPatchField, pointMesh>&
            ) const;

            //- Reconstruction of a volScalarField to surfaceScalarField
            void reconstruct
            (
                GeometricField<scalar, fvPatchField, volMesh>&,
                const volVectorField&,
                GeometricField<scalar, fvsPatchField, surfaceMesh>&,
                GeometricField<scalar, fvsPatchField, surfaceMesh>&
            );

            //- Reconstruction of a volVectorField to surfaceVectorField
            void reconstruct
            (
                volVectorField&,
                const volTensorField&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&,
                GeometricField<vector, fvsPatchField, surfaceMesh>&
            );

            //- Reconstruction for a volTensorField to surfaceTensorField
            void reconstruct
            (
                volTensorField&,
                const volTensorField&,
                const volTensorField&,
                const volTensorField&,
                GeometricField<tensor, fvsPatchField, surfaceMesh>&,
                GeometricField<tensor, fvsPatchField, surfaceMesh>&
            );

        virtual bool movePoints();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
