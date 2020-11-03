/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
18-08-2020 Jeff Heylmun:    | MUSCL reconstruction
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
    Foam::MUSCLReconstructionScheme

Description
    Base class for MUSCL reconstruction

    Example of MUSCL reconstruction specification:
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
    }
    \endverbatim


SourceFiles
    MUSCLReconstructionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef MUSCLReconstructionScheme_H
#define MUSCLReconstructionScheme_H

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
                 Class MUSCLReconstructionScheme Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MUSCLReconstructionScheme
{
protected:

    //- Const reference to mesh
    const fvMesh& mesh_;

    //- Reference to fields to interpolate
    const GeometricField<Type, fvPatchField, volMesh>& phi_;

    //- Calculate the limiter
    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    calcLimiter(const scalar& dir) const = 0;


public:

    //- Runtime type information
    TypeName("MUSCLReconstructionScheme");

    // Declare run-time constructor selection tables
        declareRunTimeSelectionTable
        (
            autoPtr,
            MUSCLReconstructionScheme,
            linearMesh,
            (
                const GeometricField<Type, fvPatchField, volMesh>& phi,
                Istream& is
            ),
            (phi, is)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            MUSCLReconstructionScheme,
            quadraticMesh,
            (
                const GeometricField<Type, fvPatchField, volMesh>& phi,
                Istream& is
            ),
            (phi, is)
        );



    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        MUSCLReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        )
        :
            mesh_(phi.mesh()),
            phi_(phi)
        {}


    // Selectors

        //- Return new tmp interpolation scheme
        static autoPtr<MUSCLReconstructionScheme<Type>> New
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            const word& name,
            const word& phaseName = word::null
        );

        //- Disallow default bitwise copy construction
        MUSCLReconstructionScheme(const MUSCLReconstructionScheme&) = delete;

    //- Destructor
    virtual ~MUSCLReconstructionScheme()
    {}

    // Member Functions

        //- Return the owner and neighbor interpolated fields
        void
        interpolateOwnNei
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


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MUSCLReconstructionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define defineMUSCLReconstructionTypes(ORDER, TYPE)                            \
typedef ORDER##MUSCLReconstructionScheme<TYPE>                                 \
    ORDER##MUSCLReconstructionScheme##TYPE;                                    \
defineNamedTemplateTypeNameAndDebug(ORDER##MUSCLReconstructionScheme##TYPE, 0);\
defineTemplateRunTimeSelectionTable                                            \
(                                                                              \
    MUSCLReconstructionScheme<TYPE>,                                           \
    ORDER##Mesh                                                                \
);


#define makeMUSCLReconstructionSchemeType(ORDER, SS, Type)                     \
                                                                               \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
                                                                               \
MUSCLReconstructionScheme<Type>::add##ORDER##MeshConstructorToTable<SS<Type>>           \
    add##ORDER##SS##Type##ORDER##MeshConstructorToTable_;                             \
                                                                               \

#define makeMUSCLReconstructionSchemes(ORDER, SS)                              \
                                                                               \
makeMUSCLReconstructionSchemeType(ORDER, SS, scalar)                           \
makeMUSCLReconstructionSchemeType(ORDER, SS, vector)                           \
makeMUSCLReconstructionSchemeType(ORDER, SS, sphericalTensor)                  \
makeMUSCLReconstructionSchemeType(ORDER, SS, symmTensor)                       \
makeMUSCLReconstructionSchemeType(ORDER, SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MUSCLReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
