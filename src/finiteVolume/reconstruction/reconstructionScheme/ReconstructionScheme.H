/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2024
     \\/     M anipulation  | Synthetik Applied Technology
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
    Foam::ReconstructionScheme

Description
    Base class for own/nei reconstruction

    Example of reconstruction specification:
    \verbatim
    interpolationSchemes
    {
        defualt linear;

        // linear MUSCL with Minmod limiter
        reconstruct(rho) linearMUSCL Minmod;

        // quadratic with Minmod limiter
        reconstruct(U)   quadraticMUSCL Minmod;

        // upwind MUSCL
        reconstruct(p)   upwindMUSCL;

        // Standard linear interpolation with Minmod limiter
        reconstruct(e)   Minmod;

        // No own and neighbour reconstruction use same field
        reconstruct(alpha)   none;
        interpolate(alpha)   CICSAM 0.5;
    }
    \endverbatim


SourceFiles
    ReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef ReconstructionScheme_H
#define ReconstructionScheme_H

#include "reconstruction.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class ReconstructionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class ReconstructionScheme
:
    public reconstruction
{
protected:

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Reference to fields to interpolate
    const GeometricField<Type, fvPatchField, volMesh>& phi_;

    //- Overwrite existing data
    bool overwrite_;

    //- Calculate the limiter
    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    calcLimiter(const scalar& dir) const = 0;


public:

    //- Runtime type information
    TypeName("ReconstructionScheme");

    // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            ReconstructionScheme,
            dictionary,
            (
                const GeometricField<Type, fvPatchField, volMesh>& phi,
                Istream& is,
                const bool overwrite
            ),
            (phi, is, overwrite)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            ReconstructionScheme,
            LinearMesh,
            (
                const GeometricField<Type, fvPatchField, volMesh>& phi,
                Istream& is,
                const bool overwrite
            ),
            (phi, is, overwrite)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            ReconstructionScheme,
            QuadraticMesh,
            (
                const GeometricField<Type, fvPatchField, volMesh>& phi,
                Istream& is,
                const bool overwrite
            ),
            (phi, is, overwrite)
        );



    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        ReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is,
            const bool overwrite
        )
        :
            mesh_(phi.mesh()),
            phi_(phi),
            overwrite_
            (
                overwrite
             || !phi.mesh().template foundObject
                <
                    GeometricField<Type, fvsPatchField, surfaceMesh>
                >(ownName())
             || !phi.mesh().template foundObject
                <
                    GeometricField<Type, fvsPatchField, surfaceMesh>
                >(neiName())
            )
        {
            if (overwrite != overwrite_)
            {
                DebugInfo
                    << "Attempt to used cached fields "
                    << ownName() << " " << neiName() << ": "
                    <<  (
                            !phi.mesh().template foundObject
                            <
                                GeometricField<Type, fvsPatchField, surfaceMesh>
                            >(ownName())
                          ? word(ownName() + " not found")
                          : word(ownName() + " found")
                        )
                    <<  ", "
                    <<  (
                            !phi.mesh().template foundObject
                            <
                                GeometricField<Type, fvsPatchField, surfaceMesh>
                            >(neiName())
                          ? word(neiName() + " not found")
                          : word(neiName() + " found")
                        ) << endl;
            }
        }


    // Selectors

        //- Return new tmp interpolation scheme
        static autoPtr<ReconstructionScheme<Type>> New
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const word& name,
            const bool overwrite
        );

        static autoPtr<ReconstructionScheme<Type>> New
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const word& name,
            const word& phaseName,
            const bool overwrite
        );

        //- Disallow default bitwise copy construction
        ReconstructionScheme(const ReconstructionScheme&) = delete;

    //- Destructor
    virtual ~ReconstructionScheme()
    {}

    // Member Functions

        using reconstruction::ownName;
        using reconstruction::neiName;

        //- Return the name of the owner field
        word ownName() const;

        //- Return the name of the neighbour field
        word neiName() const;

        //- Return if pure upwinding is used
        virtual bool upwind() const
        {
            return false;
        }

        //- Return the owner and neighbor interpolated fields
        virtual void interpolateOwnNei
        (
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&,
            tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>&
        ) const;

        //- Return the owner interpolated field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolateOwn() const = 0;

        //- Return the neighbor interpolated field
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolateNei() const = 0;

        //- Return the upwind interpolated field based on face flux
        virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
        interpolate
        (
            const surfaceScalarField& faceFlux
        ) const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ReconstructionScheme&) = delete;
};


template<class LimiterFunc>
class noneLimiter
:
    public LimiterFunc
{
public:

    noneLimiter(Istream& is)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        return 1;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define declareReconstructionScheme(SCHEME, TYPE)                              \
typedef SCHEME##ReconstructionScheme<TYPE> SCHEME##ReconstructionScheme##TYPE;


#define defineReconstructionTable(TABLE, SCHEME, TYPE)                         \
declareReconstructionScheme(SCHEME, TYPE)                                      \
defineNamedTemplateTypeNameAndDebug(SCHEME##ReconstructionScheme##TYPE, 0);    \
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    ReconstructionScheme<TYPE>,                                                \
    TABLE##Mesh                                                                \
);


#define makeReconstructionSchemeType(SCHEME, SS, Type)                         \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
ReconstructionScheme<Type>::add##SCHEME##MeshConstructorToTable<SS<Type>>      \
    add##SCHEME##SS##Type##SCHEME##MeshConstructorToTable_;                    \
ReconstructionScheme<Type>::adddictionaryConstructorToTable<SS<Type>>          \
    add##SCHEME##SS##Type##dictionaryConstructorToTable_;

#define makeReconstructionSchemes(SCHEME, SS)                                  \
                                                                               \
makeReconstructionSchemeType(SCHEME, SS, scalar)                               \
makeReconstructionSchemeType(SCHEME, SS, vector)                               \
makeReconstructionSchemeType(SCHEME, SS, sphericalTensor)                      \
makeReconstructionSchemeType(SCHEME, SS, symmTensor)                           \
makeReconstructionSchemeType(SCHEME, SS, tensor)




#define declareReconstructionSchemes(SCHEME)                                    \
                                                                               \
declareReconstructionScheme(SCHEME, scalar)                                     \
declareReconstructionScheme(SCHEME, vector)                                     \
declareReconstructionScheme(SCHEME, sphericalTensor)                            \
declareReconstructionScheme(SCHEME, symmTensor)                                 \
declareReconstructionScheme(SCHEME, tensor)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
