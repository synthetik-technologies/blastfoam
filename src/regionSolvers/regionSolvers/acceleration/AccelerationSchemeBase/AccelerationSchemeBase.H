/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2022
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::AccelerationSchemeBase

Description
    Base class for acceleration

SourceFiles
    AccelerationSchemeBase.C
    AccelerationSchemeBase.C

\*---------------------------------------------------------------------------*/

#ifndef AccelerationSchemeBase_H
#define AccelerationSchemeBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "accelerationScheme.H"
#include "PatchFieldSelector.H"
#include "runTimeSelectionTables.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AccelerationSchemeBase Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class Patch, class Mesh>
class AccelerationSchemeBase
:
    public accelerationScheme
{
protected:
// Protected data

        //- Field to relax
        GeometricField<Type, Patch, Mesh>& field_;

        //- PatchField acceleration
        autoPtr<PatchFieldSelector<Type>> selector_;

public:

    TypeName("AccelerationSchemeBase");

     // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            accelerationScheme,
            dictionary,
            (
                GeometricField<Type, Patch, Mesh>& field,
                autoPtr<PatchFieldSelector<Type>> selector,
                const dictionary& dict
            ),
            (field, selector, dict)
        );

    // Constructor
    AccelerationSchemeBase
    (
        const word& type,
        GeometricField<Type, Patch, Mesh>& field,
        autoPtr<PatchFieldSelector<Type>> selector,
        const dictionary& dict
    );


    // Selector
        static autoPtr<accelerationScheme> New
        (
            GeometricField<Type, Patch, Mesh>& field,
            autoPtr<PatchFieldSelector<Type>> selector,
            const dictionary& dict
        );


    //- Destructor
    virtual ~AccelerationSchemeBase();


    // Member functions

        //- Update residuals
        virtual void updateError();

        //- Return the patch size
        virtual label size() const
        {
            return field_.boundaryField()[this->patchi_].size();
        }

        //- Store the previous iteration for the give field
        virtual void storePrevIter()
        {
            selector_->storePrevIter();
        }
};

#define makeAccelerationScheme(Type, Patch, Mesh)                              \
    typedef AccelerationSchemeBase<Type, Patch, Mesh>                          \
        Type##Mesh##AccelerationScheme;                                        \
    defineTemplateTypeNameAndDebug                                             \
    (                                                                          \
        Type##Mesh##AccelerationScheme,                                        \
        0                                                                      \
    );                                                                         \
    defineTemplateRunTimeSelectionTable                                        \
    (                                                                          \
        Type##Mesh##AccelerationScheme,                                    \
        dictionary                                                             \
    );


#define makeAccelerationSchemes(Type)                                          \
    makeAccelerationScheme(Type, fvPatchField, volMesh);                       \
    makeAccelerationScheme(Type, fvsPatchField, surfaceMesh);                  \
    makeAccelerationScheme(Type, pointPatchField, pointMesh);


#define makeAccelerationSchemeType(AccelerationSchemeType, Type, Patch, Mesh)  \
    typedef AccelerationSchemeType<Type, Patch, Mesh>                          \
        AccelerationSchemeType##Type##Mesh##AccelerationScheme;                \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        AccelerationSchemeType##Type##Mesh##AccelerationScheme,                \
        0                                                                      \
    );                                                                         \
    AccelerationSchemeBase<Type, Patch, Mesh>::adddictionaryConstructorToTable \
    <                                                                          \
        AccelerationSchemeType##Type##Mesh##AccelerationScheme                 \
    > add##AccelerationSchemeType##Type##AccelerationScheme##Mesh##dictinoary##ConstructorToTable_;


#define makeAccelerationSchemeTypes(AccelerationSchemeType, Type)              \
    makeAccelerationSchemeType                                                 \
    (                                                                          \
        AccelerationSchemeType,                                                \
        Type,                                                                  \
        fvPatchField,                                                          \
        volMesh                                                                \
    );                                                                         \
    makeAccelerationSchemeType                                                 \
    (                                                                          \
        AccelerationSchemeType,                                                \
        Type,                                                                  \
        fvsPatchField,                                                         \
        surfaceMesh                                                            \
    );                                                                         \
    makeAccelerationSchemeType                                                 \
    (                                                                          \
        AccelerationSchemeType,                                                \
        Type,                                                                  \
        pointPatchField,                                                       \
        pointMesh                                                              \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AccelerationSchemeBase.C"
    #include "AccelerationSchemeBaseNew.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
