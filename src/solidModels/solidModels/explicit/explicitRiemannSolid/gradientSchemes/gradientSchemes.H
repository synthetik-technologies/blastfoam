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

        //- Reference to displacement field
        const volVectorField& D_;

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


        //- Return the scalar limiter
        scalar limiter
        (
            const scalar& deltaMax,
            const scalar& deltaMin,
            const scalar& delta
        ) const
        {
            if (mag(delta) < small)
            {
                return 1.0;
            }
            else if (delta > 0)
            {
                return min(deltaMax/delta, 1.0);
            }
            else
            {
                return min(deltaMin/delta, 1.0);
            }
        }


        //- Return the vector limiter
        vector limiter
        (
            const vector& deltaMax,
            const vector& deltaMin,
            const vector& delta
        ) const
        {
            return vector
            (
                limiter(deltaMax[0], deltaMin[0], delta[0]),
                limiter(deltaMax[1], deltaMin[1], delta[1]),
                limiter(deltaMax[2], deltaMin[2], delta[2])
            );
        }

        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>> calcLimiter
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const GeometricField
            <
                typename outerProduct<vector, Type>::type,
                fvPatchField,
                volMesh
            >& gradPhi
        ) const
        {
            typedef Field<Type> FieldType;
            typedef GeometricField<Type, fvPatchField, volMesh> GeoFieldType;

            typedef Field<typename outerProduct<vector, Type>::type>
                GradFieldType;

            tmp<GeoFieldType> tphiLimiter
            (
                GeoFieldType::New
                (
                    "phiLimiter",
                    mesh_,
                    dimensioned<Type>("pLimiter", dimless, pTraits<Type>::one)
                )
            );
// //             return tphiLimiter;
            GeoFieldType& phiLimiter(tphiLimiter.ref());

            // Collect min and max values from neighbourhood
            GeoFieldType phiMin("phiMin", phi);
            GeoFieldType phiMax("phiMax", phi);

            const labelUList& owner = mesh_.owner();
            const labelUList& neighbour = mesh_.neighbour();

            forAll (owner, facei)
            {
                const label own = owner[facei];
                const label nei = neighbour[facei];

                // min values
                phiMin[own] = min(phiMin[own], phi[nei]);
                phiMin[nei] = min(phiMin[nei], phi[own]);

                // max values
                phiMax[own] = max(phiMax[own], phi[nei]);
                phiMax[nei] = max(phiMax[nei], phi[own]);
            }

            // Coupled boundaries
            forAll (phi.boundaryField(), patchi)
            {
                if (phi.boundaryField()[patchi].coupled())
                {
                    const Field<Type> pNei
                    (
                        phi.boundaryField()[patchi].patchNeighbourField()
                    );

                    const labelList& fc =
                        phi.boundaryField()[patchi].patch().faceCells();

                    forAll (fc, facei)
                    {
                        const label celli = fc[facei];

                        // min value
                        phiMin[celli] = min(phiMin[celli], pNei[facei]);

                        // max value
                        phiMax[celli] = max(phiMax[celli], pNei[facei]);
                    }
                }
            }

            // Calculate limiter

            // Get geometrical information

            const volVectorField& C = mesh_.C();
            const surfaceVectorField& Cf = mesh_.Cf();

            // Compute limiter values, internal faces
            forAll (owner, facei)
            {
                const label own = owner[facei];
                const label nei = neighbour[facei];

                vector deltaOwn(Cf[facei] - C[own]);
                vector deltaNei(Cf[facei] - C[nei]);

                phiLimiter[own] =
                    min
                    (
                        phiLimiter[own],
                        limiter
                        (
                            phiMax[own] - phi[own],
                            phiMin[own] - phi[own],
                            gradPhi[own] & deltaOwn
                        )
                    );

                phiLimiter[nei] =
                    min
                    (
                        phiLimiter[nei],
                        limiter
                        (
                            phiMax[nei] - phi[nei],
                            phiMin[nei] - phi[nei],
                            gradPhi[nei] & deltaNei
                        )
                    );
            }

            // Coupled boundaries
            forAll (phi.boundaryField(), patchi)
            {
                if (phi.boundaryField()[patchi].coupled())
                {
                    // Get patch
                    const fvPatch& p = phi.boundaryField()[patchi].patch();

                    const FieldType pNei
                    (
                        phi.boundaryField()[patchi].patchNeighbourField()
                    );

                    const labelList& fc = p.faceCells();

                    const vectorField deltaOwn(p.Cf() - p.Cn());

                    // Get gradients
                    const GradFieldType gradPhiOwn
                    (
                        gradPhi.boundaryField()[patchi].patchInternalField()
                    );

                    forAll (fc, facei)
                    {
                        const label& celli = fc[facei];

                        phiLimiter[celli] =
                            min
                            (
                                phiLimiter[celli],
                                limiter
                                (
                                    phiMax[celli] - phi[celli],
                                    phiMin[celli] - phi[celli],
                                    gradPhiOwn[facei] & deltaOwn[facei]
                                )
                            );
                    }
                }
            }

            // Do parallel communication to correct limiter on
            // coupled boundaries
            phiLimiter.correctBoundaryConditions();

            return tphiLimiter;
        }

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
        gradientSchemes(const volVectorField&);


    // Destructor
    virtual ~gradientSchemes();


    // Member functions

        // Edit

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
