/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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
    Foam::fvTimeIntegrator

Description
    Time integration for a fvMesh

SourceFiles
    fvTimeIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef fvTimeIntegrator_H
#define fvTimeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrator.H"
#include "PtrList.H"
#include "HashPtrTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "UautoPtr.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "hashedWordList.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeIntegrator Declaration
\*---------------------------------------------------------------------------*/

#define FieldName(Geo, Type) \
    CAT3(Geo, CAPITALIZE(Type), Field)

#define FieldVarName(Geo, Type, OD) \
    CAT4(Geo, CAPITALIZE(Type), Fields, OD##_)

#define OldFieldVarName(Geo, Type) FieldVarName(Geo, Type, Old)

#define DeltaFieldVarName(Geo, Type) FieldVarName(Geo, Type, Delta)

class fvTimeIntegrator
:
    public timeIntegrator
{
protected:
// Protected variables

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Current step old volume
    autoPtr<volScalarField::Internal> V0Ptr_;

    //- Current step volume
    autoPtr<volScalarField::Internal> VPtr_;

    //- FvModels
    mutable UautoPtr<fvModels> modelsPtr_;

    //- fvConstraints
    mutable UautoPtr<const fvConstraints> constraintsPtr_;

    //- List of fields that have models or constraints
    mutable hashedWordList solveFields_;

    mutable wordHashSet savedOldFields_;

    //- Stored old fields
    #define decalareOldFieldTypeTable(Type, Geo)            \
    mutable HashPtrTable                                    \
    <                                                       \
        PtrList<FieldName(Geo, Type)>                       \
    > OldFieldVarName(Geo, Type);

    FOR_ALL_FIELD_TYPES(decalareOldFieldTypeTable, vol);
    FOR_ALL_FIELD_TYPES(decalareOldFieldTypeTable, surface);
    FOR_ALL_FIELD_TYPES(decalareOldFieldTypeTable, point);

    #undef decalareOldFieldTypeTable


    //- Stored delta fields
    #define decalareDeltaFieldTypeTable(Type, Geo)          \
    mutable HashPtrTable                                    \
    <                                                       \
        PtrList<FieldName(Geo, Type)>                       \
    > DeltaFieldVarName(Geo, Type);

    FOR_ALL_FIELD_TYPES(decalareDeltaFieldTypeTable, vol);
    FOR_ALL_FIELD_TYPES(decalareDeltaFieldTypeTable, surface);
    FOR_ALL_FIELD_TYPES(decalareDeltaFieldTypeTable, point);

    #undef decalareDeltaFieldTypeTable



    // Private member functions

        // Update coefficients
        virtual void update();

        //- Update all stored systems
        virtual void updateAll();

        //- Post update of all stored systems
        virtual void postUpdateAll();

        //- Insert a list of old fields into the given hash table
        template<class FieldType>
        void insertOldList
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Insert a list of delta fields into the given hash table
        template<class FieldType>
        void insertDeltaList
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Return the List of old fields
        template<class FieldType>
        const PtrList<FieldType>& lookupOld
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Return the List of old fields
        template<class FieldType>
        PtrList<FieldType>& lookupOldRef
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Return the List of delta fields
        template<class FieldType>
        const PtrList<FieldType>& lookupDelta
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Return the List of delta fields
        template<class FieldType>
        PtrList<FieldType>& lookupDeltaRef
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Clear and resize stored old fields
        template<class FieldType>
        void clearOldFields(HashPtrTable<PtrList<FieldType>>& table);

        template<class FieldType>
        void clearDeltaFields(HashPtrTable<PtrList<FieldType>>& table);

        //- Reset all fields of the given type to the oldTime
        template<class FieldType>
        void resetFields();

        //- Reset all fields of the given type to the oldTime
        template<class FieldType>
        void shuffleFields
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const labelList& oldToNew,
            const word& type
        );


public:

    //- Runtime type information
    TypeName("fvTimeIntegrator");

    // Constructor
    fvTimeIntegrator(const fvMesh& mesh);


    //- Destructor
    virtual ~fvTimeIntegrator();


    // Selectors

        //- Return a reference to a timeIntegrator object
        static fvTimeIntegrator& NewRef(const fvMesh& mesh);


    // Member Functions

        //- Update anything that needs to be done before mesh updates
        void preUpdateMesh();

        //- Integrate fluxes in time
        virtual void integrate();

        //- Clear all ODE fields (old/delta)
        virtual void clear();

        //- Reset fields to the old time
        virtual void reset();

        //- Create/lookup fvModels/fvConstraints
        void createModels() const;

        //- Return the fvModels
        inline const fvModels& models() const
        {
            if (!modelsPtr_.valid())
            {
                createModels();
            }
            return modelsPtr_();
        }

        //- Access the fvModels
        inline fvModels& modelsRef() const
        {
            if (!modelsPtr_.valid())
            {
                createModels();
            }
            return modelsPtr_();
        }

        //- Return the fvConstraints
        inline const fvConstraints& constraints() const
        {
            if (!constraintsPtr_.valid())
            {
                createModels();
            }
            return constraintsPtr_();
        }

        //- Is a field modified by fvModels or fvConstraints
        inline bool needSolve(const word& field) const
        {
            if (!modelsPtr_.valid())
            {
                createModels();
            }
            return solveFields_.found(field);
        }

        //- Return current step old volume
        virtual tmp<scalarField> V0() const;

        //- Return current step volume
        virtual tmp<scalarField> V() const;

        //- Total old volume
        virtual scalar totalV0() const;

        //- Total volume
        virtual scalar totalV() const;

        //- Add fields to saved fields
        #define decalareAddFieldType(Type, Geo)                 \
        inline void addOldField                                 \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;                                                \
        inline void addDeltaField                               \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;
        FOR_ALL_FIELD_TYPES(decalareAddFieldType, vol);
        FOR_ALL_FIELD_TYPES(decalareAddFieldType, surface);
        FOR_ALL_FIELD_TYPES(decalareAddFieldType, point);

        #undef decalareAddFieldType

        //- Lookup fields
        #define decalareLookupFieldType(Type, Geo)              \
        inline const PtrList<FieldName(Geo, Type)>&             \
        oldFields                                               \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;                                                \
                                                                \
        inline PtrList<FieldName(Geo, Type)>&                   \
        oldFieldsRef                                            \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;                                                \
                                                                \
        inline const PtrList<FieldName(Geo, Type)>&             \
        deltaFields                                             \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;                                                \
                                                                \
        inline PtrList<FieldName(Geo, Type)>&                   \
        deltaFieldsRef                                          \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;
        FOR_ALL_FIELD_TYPES(decalareLookupFieldType, vol);
        FOR_ALL_FIELD_TYPES(decalareLookupFieldType, surface);
        FOR_ALL_FIELD_TYPES(decalareLookupFieldType, point);

        #undef decalareLookupFieldType


        //- Dummy IO operations
        virtual bool writeData(Ostream& os) const
        {
            return true;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvTimeIntegratorI.H"

#ifdef NoRepository
#include "fvTimeIntegratorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
