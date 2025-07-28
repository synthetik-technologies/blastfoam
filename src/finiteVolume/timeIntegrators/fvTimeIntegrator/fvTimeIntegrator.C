/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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

\*---------------------------------------------------------------------------*/

#include "fvTimeIntegrator.H"
#include "timeIntegrationSystem.H"
#include "pointFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvTimeIntegrator, 0);
}


// * * * * * * * * * * * * * * Protected Functions * * * * * * * * * * * * * //

void Foam::fvTimeIntegrator::update()
{
    if (coeffs_->update(curTimeIndex_, time().value()))
    {
        initialize();
    }
    const labelList oldMap(coeffs_->oldSaveMap());
    if (oldMap.size())
    {
        #define ShuffleOldFieldTypes(Type, map, type)                      \
            shuffleFields<FieldName(vol, Type)>                            \
            (                                                              \
                OldFieldVarName(vol, Type),                                \
                map,                                                       \
                type                                                       \
            );                                                             \
            shuffleFields<FieldName(surface, Type)>                        \
            (                                                              \
                OldFieldVarName(surface, Type),                            \
                map,                                                       \
                type                                                       \
            );                                                             \
            shuffleFields<FieldName(point, Type)>                          \
            (                                                              \
                OldFieldVarName(point, Type),                              \
                map,                                                       \
                type                                                       \
            );
        FOR_ALL_FIELD_TYPES(ShuffleOldFieldTypes, oldMap, "old");
        #undef ShuffleOldFieldTypes
    }

    const labelList deltaMap(coeffs_->deltaSaveMap());
    if (deltaMap.size())
    {
        #define ShuffleDeltaFieldTypes(Type, map, type)                    \
            shuffleFields<FieldName(vol, Type)>                            \
            (                                                              \
                DeltaFieldVarName(vol, Type),                              \
                map,                                                       \
                type                                                       \
            );                                                             \
            shuffleFields<FieldName(surface, Type)>                        \
            (                                                              \
                DeltaFieldVarName(surface, Type),                          \
                map,                                                       \
                type                                                       \
            );                                                             \
            shuffleFields<FieldName(point, Type)>                          \
            (                                                              \
                DeltaFieldVarName(point, Type),                            \
                map,                                                       \
                type                                                       \
            );
        FOR_ALL_FIELD_TYPES(ShuffleDeltaFieldTypes, deltaMap, "delta");
        #undef ShuffleDeltaFieldTypes
    }

    forAll(systems_, i)
    {
        systems_[i].save();
    }
}


void Foam::fvTimeIntegrator::updateAll()
{
    // Set use a linear change in volume
    // All fields are scaled according to the true volume
    if (mesh_.moving())
    {
        tmp<volScalarField::Internal> tV0(mesh_.V0());
        const volScalarField::Internal& V0 = tV0();

        tmp<volScalarField::Internal> tV(mesh_.V());
        const volScalarField::Internal& V = tV();

        if (!V0ByVPtr_.valid())
        {
            V0ByVPtr_.set
            (
                new volScalarField::Internal
                (
                    IOobject
                    (
                        "fvTimeIntegrator:V0ByV",
                        mesh_.time().name(),
                        mesh_
                    ),
                    V0/V
                )
            );
        }
        else
        {
            V0ByVPtr_() = V0/V;
        }


//         if (!VPtr_.valid())
//         {
//             VPtr_.set
//             (
//                 new volScalarField::Internal
//                 (
//                     IOobject
//                     (
//                         "fvTimeIntegrator:V",
//                         mesh_.time().name(),
//                         mesh_
//                     ),
//                     mesh_,
//                     1.0//V0 + f()*(V - V0)
//                 )
//             );
//         }
//         else
//         {
//             VPtr_() = 1.0;//V0 + f()*(V - V0);
//         }

//         if (!V0Ptr_.valid())
//         {
//             V0Ptr_.set
//             (
//                 new volScalarField::Internal
//                 (
//                     IOobject
//                     (
//                         "fvTimeIntegrator:V0",
//                         mesh_.time().name(),
//                         mesh_
//                     ),
//                     V0ByVPtr_()()*VPtr_()
// //                     (V0 + f0()*(V - V0))
//                 )
//             );
//         }
//         else
//         {
//             V0Ptr_() = V0ByVPtr_()()*VPtr_();//(V0 + f0()*(V - V0));
//         }
    }
    forAll(systems_, i)
    {
        systems_[i].update();
    }
}


void Foam::fvTimeIntegrator::postUpdateAll()
{
    forAll(systems_, i)
    {
        Info<< "Post-updating " << systems_[i].name() << ":" << endl;
        systems_[i].postUpdate();
    }
    Info<< endl;

    if (modelsPtr_.valid())
    {
        modelsPtr_->correct();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvTimeIntegrator::fvTimeIntegrator(const fvMesh& mesh)
:
    timeIntegrator(mesh, mesh.schemes().dict().subDict("ddtSchemes")),
    mesh_(mesh),
    V0Ptr_(nullptr),
    VPtr_(nullptr),
    V0ByVPtr_(nullptr),
    modelsPtr_(nullptr),
    constraintsPtr_(nullptr),
    solveFields_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvTimeIntegrator::~fvTimeIntegrator()
{}


// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //


void Foam::fvTimeIntegrator::createModels() const
{
    DebugInfo<< "Creating fvModels and fvConstraints" << endl;
    modelsPtr_.set(&fvModels::New(const_cast<fvMesh&>(mesh_)));
    constraintsPtr_.set(&fvConstraints::New(mesh_));

    const PtrList<fvModel>& models(modelsPtr_());
    forAll(models, modeli)
    {
        wordList fields(models[modeli].addSupFields());
        forAll(fields, fieldi)
        {
            if (!solveFields_.found(fields[fieldi]))
            {
                solveFields_.append(fields[fieldi]);
            }
        }
    }

    const PtrList<fvConstraint>& constraints(constraintsPtr_());
    forAll(constraints, modeli)
    {
        wordList fields(constraints[modeli].constrainedFields());
        forAll(fields, fieldi)
        {
            if (!solveFields_.found(fields[fieldi]))
            {
                solveFields_.append(fields[fieldi]);
            }
        }
    }
    DebugInfo<< "Fields to solve:" << nl<< solveFields_ << endl;
}


void Foam::fvTimeIntegrator::preUpdateMesh()
{
    DebugInfo<< "Post Update" << endl;
    forAll(systems_, i)
    {
        systems_[i].preUpdateMesh();
    }

    if (modelsPtr_.valid())
    {
        modelsPtr_->preUpdateMesh();
    }
}


void Foam::fvTimeIntegrator::integrate()
{
    timeIntegrator::integrate();

    if (!obr_.time().subCycling())
    {
        DebugInfo<< "Clearing SourceTerms fields" << endl;
        #define ClearSourceTypes(Type, Geo)             \
            FieldVarName(Geo, Type, Source).clear();    \
            FieldVarName(Geo, Type, IntegratedSource).clear();
        FOR_ALL_FIELD_TYPES(ClearSourceTypes, vol);
        #undef ClearSourceTypes
    }
}


void Foam::fvTimeIntegrator::clear()
{
    if (!coeffs_->save())
    {
        DebugInfo<< "Clearing ODE fields" << endl;
        forAll(systems_, i)
        {
            systems_[i].clear();
        }

        #define ClearFieldTypes(Type, Geo)                  \
            clearOldFields<FieldName(Geo, Type)>            \
            (                                               \
                OldFieldVarName(Geo, Type)                  \
            );                                              \
            clearDeltaFields<FieldName(Geo, Type)>          \
            (                                               \
                DeltaFieldVarName(Geo, Type)                \
            );
        FOR_ALL_FIELD_TYPES(ClearFieldTypes, vol);
        FOR_ALL_FIELD_TYPES(ClearFieldTypes, surface);
        FOR_ALL_FIELD_TYPES(ClearFieldTypes, point);
        #undef ClearFieldTypes
    }

    DebugInfo<< "Clearing SourceTerms fields" << endl;
    #define ClearSourceTypes(Type, Geo)             \
        FieldVarName(Geo, Type, Source).clear();    \
        FieldVarName(Geo, Type, IntegratedSource).clear();
    FOR_ALL_FIELD_TYPES(ClearSourceTypes, vol);
    #undef ClearSourceTypes
}


void Foam::fvTimeIntegrator::reset()
{
    DebugInfo<< "Resetting fields to old time" << endl;
    #define ResetOldFieldTypes(Type, Geo)               \
        resetFields<FieldName(Geo, Type)>();
    FOR_ALL_FIELD_TYPES(ResetOldFieldTypes, vol);
    FOR_ALL_FIELD_TYPES(ResetOldFieldTypes, surface);
    FOR_ALL_FIELD_TYPES(ResetOldFieldTypes, point);
    #undef ResetOldFieldTypes
}


Foam::tmp<Foam::scalarField> Foam::fvTimeIntegrator::V0() const
{
    if (stepi_ == 0 && mesh_.moving())
    {
        return V0ByVPtr_();
        // return mesh_.V0().primitiveField()/mesh_.V().primitiveField();
    }
    return tmp<scalarField>();
}


Foam::tmp<Foam::scalarField> Foam::fvTimeIntegrator::V() const
{
//     if (VPtr_.valid())
//     {
//         return VPtr_();
//     }
    return tmp<scalarField>();
}


Foam::scalar Foam::fvTimeIntegrator::totalV0() const
{
    if (stepi_ == 0 && mesh_.moving())
    {
        return gSum(mesh_.V0());
    }
    return gSum(mesh_.V());
}

Foam::scalar Foam::fvTimeIntegrator::totalV() const
{
    return gSum(mesh_.V());
}

#define defineSourceLookupType(Type, Geo)                                   \
bool Foam::fvTimeIntegrator::foundSource                                    \
(                                                                           \
    const FieldName(Geo, Type)& f                                           \
) const                                                                     \
{                                                                           \
    return                                                                  \
        FieldVarName(Geo, Type, IntegratedSource).found(f.name())           \
     || FieldVarName(Geo, Type, Source).found(f.name());                    \
}                                                                           \
                                                                            \
void Foam::fvTimeIntegrator::addSource                                      \
(                                                                           \
    const word& fName,                                                      \
    const FieldName(Geo, Type)::Internal& S                                 \
)                                                                           \
{                                                                           \
    if (!mesh_.foundObject<FieldName(Geo, Type)::Internal>(fName))          \
    {                                                                       \
        FatalErrorInFunction                                                \
            << fName  << " is not a registered "                            \
            << FieldName(Geo, Type)::Internal::typeName << endl             \
            << abort(FatalError);                                           \
    }                                                                       \
    HashPtrTable<FieldName(Geo, Type)::Internal>::iterator iter =           \
        FieldVarName(Geo, Type, Source).find(fName);                        \
    if (iter != FieldVarName(Geo, Type, Source).end())                      \
    {                                                                       \
        (*iter()) += S;                                                     \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        FieldVarName(Geo, Type, Source).insert                              \
        (                                                                   \
            fName,                                                          \
            new FieldName(Geo, Type)::Internal(fName + "Source" , S)        \
        );                                                                  \
    }                                                                       \
}                                                                           \
                                                                            \
void Foam::fvTimeIntegrator::addSource                                      \
(                                                                           \
    const word& fName,                                                      \
    const tmp<FieldName(Geo, Type)::Internal>& S                            \
)                                                                           \
{                                                                           \
    if (!mesh_.foundObject<FieldName(Geo, Type)::Internal>(fName))          \
    {                                                                       \
        FatalErrorInFunction                                                \
            << fName  << " is not a registered "                            \
            << FieldName(Geo, Type)::Internal::typeName << endl             \
            << abort(FatalError);                                           \
    }                                                                       \
    HashPtrTable<FieldName(Geo, Type)::Internal>::iterator iter =           \
        FieldVarName(Geo, Type, Source).find(fName);                        \
    if (iter != FieldVarName(Geo, Type, Source).end())                      \
    {                                                                       \
        (*iter()) += S;                                                     \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        S.ref().rename(fName + "Source");                                   \
        FieldVarName(Geo, Type, Source).insert                              \
        (                                                                   \
            fName,                                                          \
            S.ptr()                                                         \
        );                                                                  \
    }                                                                       \
}                                                                           \
                                                                            \
void Foam::fvTimeIntegrator::addIntegratedSource                            \
(                                                                           \
    const word& fName,                                                      \
    const FieldName(Geo, Type)::Internal& S                                 \
)                                                                           \
{                                                                           \
    if (!mesh_.foundObject<FieldName(Geo, Type)::Internal>(fName))          \
    {                                                                       \
        FatalErrorInFunction                                                \
            << fName  << " is not a registered "                            \
            << FieldName(Geo, Type)::Internal::typeName << endl             \
            << abort(FatalError);                                           \
    }                                                                       \
    HashPtrTable<FieldName(Geo, Type)::Internal>::iterator iter =           \
        FieldVarName(Geo, Type, IntegratedSource).find(fName);              \
    if (iter != FieldVarName(Geo, Type, IntegratedSource).end())            \
    {                                                                       \
        (*iter()) += S;                                                     \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        FieldVarName(Geo, Type, IntegratedSource).insert                    \
        (                                                                   \
            fName,                                                          \
            new FieldName(Geo, Type)::Internal                              \
            (                                                               \
                fName + "IntegratedSource" ,                                \
                S                                                           \
            )                                                               \
        );                                                                  \
    }                                                                       \
}                                                                           \
                                                                            \
void Foam::fvTimeIntegrator::addIntegratedSource                            \
(                                                                           \
    const word& fName,                                                      \
    const tmp<FieldName(Geo, Type)::Internal>& S                            \
)                                                                           \
{                                                                           \
    if (!mesh_.foundObject<FieldName(Geo, Type)::Internal>(fName))          \
    {                                                                       \
        FatalErrorInFunction                                                \
            << fName  << " is not a registered "                            \
            << FieldName(Geo, Type)::Internal::typeName << endl             \
            << abort(FatalError);                                           \
    }                                                                       \
    HashPtrTable<FieldName(Geo, Type)::Internal>::iterator iter =           \
        FieldVarName(Geo, Type, IntegratedSource).find(fName);              \
    if (iter != FieldVarName(Geo, Type, IntegratedSource).end())            \
    {                                                                       \
        (*iter()) += S;                                                     \
    }                                                                       \
    else                                                                    \
    {                                                                       \
        S.ref().rename(fName + "IntegratedSource");                         \
        FieldVarName(Geo, Type, IntegratedSource).insert                    \
        (                                                                   \
            fName,                                                          \
            S.ptr()                                                         \
        );                                                                  \
    }                                                                       \
}                                                                           \
                                                                            \
bool Foam::fvTimeIntegrator::addDeltaSource                                 \
(                                                                           \
    const word& fName,                                                      \
    FieldName(Geo, Type)::Internal& fDelta                                  \
) const                                                                     \
{                                                                           \
    bool set = false;                                                       \
    {                                                                       \
        HashPtrTable<FieldName(Geo, Type)::Internal>::const_iterator iter = \
            FieldVarName(Geo, Type, Source).find(fName);                    \
        if (iter != FieldVarName(Geo, Type, Source).cend())                 \
        {                                                                   \
            fDelta += *iter();                                              \
            set = true;                                                     \
        }                                                                   \
    }                                                                       \
    {                                                                       \
        HashPtrTable<FieldName(Geo, Type)::Internal>::const_iterator iter = \
            FieldVarName(Geo, Type, IntegratedSource).find(fName);          \
        if (iter != FieldVarName(Geo, Type, IntegratedSource).cend())       \
        {                                                                   \
            const TimeState& ts =                                           \
                time().subCycling()                                         \
              ? time().prevTimeState()                                      \
              : time();                                                     \
            fDelta += (*iter())/mesh_.V()/ts.deltaT();                      \
            set = true;                                                     \
        }                                                                   \
    }                                                                       \
    return set;                                                             \
}

FOR_ALL_FIELD_TYPES(defineSourceLookupType, vol);

#undef defineSourceLookupType

// ************************************************************************* //
