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
    Foam::MUSCLReconstruction

Description
    Base class to hold limiter and calculate the limiter fields. A limiter is
    computed for each component of a field

SourceFiles
    MUSCLReconstruction.C

\*---------------------------------------------------------------------------*/

#ifndef MUSCLReconstruction_H
#define MUSCLReconstruction_H

#include "MUSCLReconstructionScheme.H"
#include "LimitFuncs.H"
#include "NVDTVD.H"
#include "NVDVTVDV.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class MUSCLReconstruction Declaration
\*---------------------------------------------------------------------------*/

template
<
    class Type,
    class MUSCLType,
    class Limiter,
    template<class> class LimitFunc
>
class MUSCLReconstruction
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
    TypeName("MUSCLReconstruction");


    // Constructors

        //- Construct from mesh and faceFlux and limiter scheme
        MUSCLReconstruction
        (
            const GeometricField<Type, fvPatchField, volMesh>& phi,
            Istream& is
        )
        :
            MUSCLType(phi, is),
            Limiter(is)
        {}

    //- Destructor
    virtual ~MUSCLReconstruction()
    {}

        //- Disallow default bitwise copy construction
        MUSCLReconstruction(const MUSCLReconstruction&) = delete;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const MUSCLReconstruction&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables
#define makeMUSCLReconstructionType\
(                                                                              \
    ORDER,                                                                     \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef MUSCLReconstruction                                                    \
    <                                                                          \
        TYPE,                                                                  \
        ORDER##MUSCLReconstructionScheme##TYPE,                                \
        LIMITER<NVDTVD>,                                                       \
        limitFuncs::LIMFUNC                                                    \
    > ORDER##MUSCLReconstruction##TYPE##LIMITER##NVDTVD##LIMFUNC##_;           \
                                                                               \
defineTemplateTypeNameAndDebugWithName                                         \
    (ORDER##MUSCLReconstruction##TYPE##LIMITER##NVDTVD##LIMFUNC##_, #SS, 0);   \
                                                                               \
MUSCLReconstructionScheme<TYPE>::add##ORDER##MeshConstructorToTable            \
    <                                                                          \
        MUSCLReconstruction                                                    \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LIMITER<NVDTVD>,                                                   \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
    > add##SS##LIMFUNC##TYPE##ORDER##MeshConstructorToTable_;


#define makeMUSCLReconstruction(ORDER, SS, LIMITER)\
                                                                               \
makeMUSCLReconstructionType(ORDER, SS, LIMITER, NVDTVD, magSqr, scalar)        \
makeMUSCLReconstructionType(ORDER, SS, LIMITER, NVDTVD, magSqr, vector)        \
makeMUSCLReconstructionType                                                    \
(                                                                              \
    ORDER,                                                                     \
    SS,                                                                        \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    magSqr,                                                                    \
    sphericalTensor                                                            \
)                                                                              \
makeMUSCLReconstructionType(ORDER, SS, LIMITER, NVDTVD, magSqr, symmTensor)    \
makeMUSCLReconstructionType(ORDER, SS, LIMITER, NVDTVD, magSqr, tensor)


#define makeLMUSCLReconstruction\
(                                                                              \
    ORDER,                                                                     \
    SS,                                                                        \
    LLIMITER,                                                                  \
    LIMITER,                                                                   \
    NVDTVD,                                                                    \
    LIMFUNC,                                                                   \
    TYPE                                                                       \
)                                                                              \
                                                                               \
typedef MUSCLReconstruction                                                    \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LLIMITER<LIMITER<NVDTVD>>,                                         \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
        ORDER##MUSCLReconstruction##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_;\
defineTemplateTypeNameAndDebugWithName                                         \
    (                                                                          \
        ORDER##MUSCLReconstruction##TYPE##LLIMITER##LIMITER##NVDTVD##LIMFUNC##_,\
        #SS,                                                                   \
        0                                                                      \
    );                                                                         \
                                                                               \
MUSCLReconstructionScheme<TYPE>::add##ORDER##MeshConstructorToTable            \
    <                                                                          \
        MUSCLReconstruction                                                    \
        <                                                                      \
            TYPE,                                                              \
            ORDER##MUSCLReconstructionScheme##TYPE,                            \
            LLIMITER<LIMITER<NVDTVD>>,                                         \
            limitFuncs::LIMFUNC                                                \
        >                                                                      \
    > add##SS##LIMFUNC##TYPE##ORDER##MeshConstructorToTable_;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MUSCLReconstruction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
