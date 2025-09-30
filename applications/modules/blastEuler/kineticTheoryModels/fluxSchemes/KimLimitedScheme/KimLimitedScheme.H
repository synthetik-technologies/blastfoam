/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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
    Foam::KimLimitedScheme

Description
    Class to create NVD/TVD limited weighting-factors.

    The particular differencing scheme class is supplied as a template
    argument, the weight function of which is called by the weight function
    of this class for the internal faces as well as faces of coupled
    patches (e.g. processor-processor patches). The weight function is
    supplied the central-differencing weighting factor, the face-flux, the
    cell and face gradients (from which the normalised variable
    distribution may be created) and the cell centre distance.

    This code organisation is both neat and efficient, allowing for
    convenient implementation of new schemes to run on parallelised cases.

SourceFiles
    KimLimitedScheme.C

\*---------------------------------------------------------------------------*/

#ifndef KimLimitedScheme_H
#define KimLimitedScheme_H

#include "limitedSurfaceInterpolationScheme.H"
#include "LimitFuncs.H"
#include "NVDTVD.H"
#include "NVDVTVDV.H"
#include "kineticTheorySystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class KimScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type, class Limiter, template<class> class LimitFunc>
class KimLimitedScheme
:
    public limitedSurfaceInterpolationScheme<Type>,
    public Limiter
{
    // Private Member data

        //- Smoothing coefficient
        scalar D_;

        //- Reference to kinetic theory system
        const kineticTheorySystem& kineticTheory_;

    // Private Member Functions

        //- Calculate the limiter
        void calcLimiter
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            surfaceScalarField& limiterField
        ) const;

        scalar limiter
        (
            const scalar cdWeight,
            const scalar faceFlux,
            const typename Limiter::phiType& phiP,
            const typename Limiter::phiType& phiN,
            const typename Limiter::gradPhiType& gradcP,
            const typename Limiter::gradPhiType& gradcN,
            const vector& d,
            const scalar G
        ) const
        {
            scalar r = Limiter::r
            (
                faceFlux, phiP, phiN, gradcP, gradcN, d
            );

            return max(0, G*min(1.0, r));
        }


public:

    //- Runtime type information
    TypeName("Kim");

    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        KimLimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            const Limiter& weight
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(weight),
            D_(1.0),
            kineticTheory_
            (
                mesh.lookupObject<kineticTheorySystem>("kineticTheorySystem")
            )
        {}

        //- Construct from mesh and Istream.
        //  The name of the flux field is read from the Istream and looked-up
        //  from the mesh objectRegistry
        KimLimitedScheme
        (
            const fvMesh& mesh,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, is),
            Limiter(),
            D_(readScalar(is)),
            kineticTheory_
            (
                mesh.lookupObject<kineticTheorySystem>("kineticTheorySystem")
            )
        {}

        //- Construct from mesh, faceFlux and Istream
        KimLimitedScheme
        (
            const fvMesh& mesh,
            const surfaceScalarField& faceFlux,
            Istream& is
        )
        :
            limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
            Limiter(),
            D_(readScalar(is)),
            kineticTheory_
            (
                mesh.lookupObject<kineticTheorySystem>("kineticTheorySystem")
            )
        {}

        //- Disallow default bitwise copy construction
        KimLimitedScheme(const KimLimitedScheme&) = delete;


    // Member Functions

        //- Return the interpolation weighting factors
        virtual tmp<surfaceScalarField> limiter
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const KimLimitedScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeKimLimitedSurfaceInterpolationTypeScheme\
(                                                                              \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef KimLimitedScheme<TYPE, NVDTVD, limitFuncs::LIMFUNC>              \
    LimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_;                          \
defineTemplateTypeNameAndDebugWithName                                         \
    (LimitedScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);                \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshConstructorToTable                    \
<KimLimitedScheme<TYPE, NVDTVD, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshConstructorToTable_;                           \
                                                                               \
surfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable                \
<KimLimitedScheme<TYPE, NVDTVD, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToTable_;                       \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshConstructorToTable             \
<KimLimitedScheme<TYPE, NVDTVD, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshConstructorToLimitedTable_;                    \
                                                                               \
limitedSurfaceInterpolationScheme<TYPE>::addMeshFluxConstructorToTable         \
<KimLimitedScheme<TYPE, NVDTVD, limitFuncs::LIMFUNC>>                    \
    add##SS##LIMFUNC##TYPE##MeshFluxConstructorToLimitedTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KimLimitedSchemeTmp.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
