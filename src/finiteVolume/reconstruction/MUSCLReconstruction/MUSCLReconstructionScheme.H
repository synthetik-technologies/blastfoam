/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
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
    Foam::MUSCLReconstructionScheme

Description
    Base class to hold limiter and calculate the limiter fields. A limiter is
    computed for each component of a field

SourceFiles
    MUSCLReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef MUSCLReconstructionScheme_H
#define MUSCLReconstructionScheme_H

#include "ReconstructionScheme.H"
#include "LimitFuncs.H"
#include "NVDTVD.H"
#include "NVDVTVDV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MUSCLReconstructionScheme Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Type,
    class MUSCLType,
    class Limiter,
    template<class> class LimitFunc
>
class MUSCLReconstructionScheme
:
    public MUSCLType,
    public Limiter
{
protected:

    //- Calculate the limiter
    virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>
    calcLimiter(const scalar& dir) const;


public:

    //- Runtime type information
    TypeName(MUSCLType::typeName + "<" + Limiter::typeName + ">");


    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        MUSCLReconstructionScheme
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        )
        :
            MUSCLType(phi, is),
            Limiter(is)
        {}

    //- Destructor
    virtual ~MUSCLReconstructionScheme()
    {}

        //- Disallow default bitwise copy construction
        MUSCLReconstructionScheme(const MUSCLReconstructionScheme&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MUSCLReconstructionScheme&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables
#define makeMUSCLReconstructionType\
(                                                                              \
    ORDER,                                                                     \
    NAME,                                                                      \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef MUSCLReconstructionScheme                                              \
    <                                                                          \
        TYPE,                                                                  \
        ORDER##MUSCLReconstructionScheme##TYPE,                                \
        LIMITER<NVDTVD>,                                                       \
        limitFuncs::LIMFUNC                                                    \
    > ORDER##MUSCLReconstructionScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_;     \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
    (                                                                          \
        ORDER##MUSCLReconstructionScheme##TYPE##LIMITER##NVDTVD##LIMFUNC##_,   \
        NAME "<" #SS ">",                                                      \
        0                                                                      \
    );                                                                         \
                                                                               \
ReconstructionScheme<TYPE>::add##ORDER##MeshConstructorToTable                 \
    <                                                                          \
        MUSCLReconstructionScheme                                              \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LIMITER<NVDTVD>,                                                   \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
    > add##SS##LIMFUNC##TYPE##ORDER##MeshConstructorToTable_;                  \
ReconstructionScheme<TYPE>::adddictionaryConstructorToTable                    \
    <                                                                          \
        MUSCLReconstructionScheme                                              \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LIMITER<NVDTVD>,                                                   \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
    > add##ORDER##SS##LIMFUNC##TYPE##dictionaryConstructorToTable_;


#define makeMUSCLReconstruction(ORDER, NAME, SS, LIMITER)                      \
                                                                               \
makeMUSCLReconstructionType(ORDER, NAME, SS, LIMITER, NVDTVD, magSqr, scalar)  \
makeMUSCLReconstructionType(ORDER, NAME, SS, LIMITER, NVDTVD, magSqr, vector)  \
makeMUSCLReconstructionType                                                    \
(                                                                              \
    ORDER,                                                                     \
    NAME,                                                                      \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    magSqr,                                                                    \
    sphericalTensor                                                            \
)                                                                              \
makeMUSCLReconstructionType(ORDER, NAME, SS, LIMITER, NVDTVD, magSqr, symmTensor) \
makeMUSCLReconstructionType(ORDER, NAME, SS, LIMITER, NVDTVD, magSqr, tensor)


#define makeLMUSCLReconstruction\
(                                                                              \
    ORDER,                                                                     \
    NAME,                                                                      \
    SS,                                                                        \
    LLIMITER,                                                                  \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef MUSCLReconstructionScheme                                                    \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LLIMITER<LIMITER<NVDTVD>>,                                         \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
        ORDER##MUSCLReconstructionScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_;\
defineTemplateTypeNameAndDebugWithName                                         \
    (                                                                          \
        ORDER##MUSCLReconstructionScheme##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_,\
        NAME "<" #SS ">",                                                      \
        0                                                                      \
    );                                                                         \
                                                                               \
ReconstructionScheme<TYPE>::add##ORDER##MeshConstructorToTable                 \
    <                                                                          \
        MUSCLReconstructionScheme                                              \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LLIMITER<LIMITER<NVDTVD>>,                                         \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
    > add##SS##LIMFUNC##TYPE##ORDER##MeshConstructorToTable_;                  \
ReconstructionScheme<TYPE>::adddictionaryConstructorToTable                    \
    <                                                                          \
        MUSCLReconstructionScheme                                              \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LLIMITER<LIMITER<NVDTVD>>,                                         \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
    > add##ORDER##SS##LIMFUNC##TYPE##dictinoaryConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MUSCLReconstructionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
