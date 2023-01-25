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
    Foam::timeIntegrator

Description
    Base class for time integration

SourceFiles
    timeIntegrator.C
    newTimeIntegrator.C

\*---------------------------------------------------------------------------*/

#ifndef timeIntegrator_H
#define timeIntegrator_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "runTimeSelectionTables.H"
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

class timeIntegrationSystem;

/*---------------------------------------------------------------------------*\
                           Class timeIntegrator Declaration
\*---------------------------------------------------------------------------*/

#define FieldName(Geo, Type) \
    CAT3(Geo, CAPITALIZE(Type), Field)

#define FieldVarName(Geo, Type, OD) \
    CAT4(Geo, CAPITALIZE(Type), Fields, OD##_)

#define OldFieldVarName(Geo, Type) FieldVarName(Geo, Type, Old)

#define DeltaFieldVarName(Geo, Type) FieldVarName(Geo, Type, Delta)

class timeIntegrator
:
    public regIOobject
{
protected:
// Protected variables

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to compressible system
    UPtrList<timeIntegrationSystem> systems_;

    //- Current step index
    label stepi_;

    //- Old value coefficients
    List<List<scalar>> as_;

    //- Delta coefficients
    List<List<scalar>> bs_;

    //- Time step fractions
    scalarList f_;

    //- Starting time step fractions
    scalarList f0_;

    //- Number of steps
    label nSteps_;

    //- Stored old field indexes
    labelList oldIs_;

    //- Number of stored fields
    label nOld_;

    //- Stored delta indexes
    labelList deltaIs_;

    //- Number of stored deltas
    label nDelta_;

    //- Current step V0/V
    volScalarField::Internal V0byV_;

    //- Current time index
    label curTimeIndex_;

    //- Is this a restarted time step
    bool restart_;

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

        //- Update all stored systems
        void updateAll();

        //- Post update of all stored systems
        void postUpdateAll();

        //- Initialize ODE sizes
        void initialize();

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
        PtrList<FieldType>& lookupOld
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        );

        //- Return the List of delta fields
        template<class FieldType>
        const PtrList<FieldType>& lookupDelta
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        ) const;

        //- Return the List of delta fields
        template<class FieldType>
        PtrList<FieldType>& lookupDelta
        (
            HashPtrTable<PtrList<FieldType>>& table,
            const FieldType& f
        );

        //- Clear and resize stored old fields
        template<class FieldType>
        void clearOldFields(HashPtrTable<PtrList<FieldType>>& table);

        template<class FieldType>
        void clearDeltaFields(HashPtrTable<PtrList<FieldType>>& table);

        //- Reset all fields of the given type to the oldTime
        template<class FieldType>
        void resetFields();

public:

    //- Runtime type information
    TypeName("timeIntegrator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            timeIntegrator,
            dictionary,
            (const fvMesh& mesh, Istream& is),
            (mesh, is)
        );

    // Constructor
    timeIntegrator(const fvMesh& mesh);


    //- Destructor
    virtual ~timeIntegrator();


    // Selectors

        //- Return a new timeIntegrator object
        static autoPtr<timeIntegrator> New(const fvMesh& mesh);

        //- Return a reference to a timeIntegrator object
        static timeIntegrator& NewRef(const fvMesh& mesh);


    // Member Functions

        //- Add a system to the
        void addSystem(timeIntegrationSystem& system);

        //- Update anything that needs to be done before mesh updates
        void preUpdateMesh();

        //- Integrate fluxes in time
        virtual void integrate();

        //- Clear all ODE fields (old/delta)
        void clear();

        //- Reset fields to the old time
        void reset();

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
        inline fvModels& models()
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

        //- Return if on a restarted time step
        bool restart() const
        {
            return restart_;
        }

        //- Return old coefficients for the current step
        inline scalarList a() const
        {
            return as_[stepi_ - 1];
        }

        //- Return delta coefficients for the current step
        inline scalarList b() const
        {
            return bs_[stepi_ - 1];
        }


        //- Return the number of steps
        inline label nSteps() const
        {
            return as_.size();
        }

        //- Return old field indexes
        inline const labelList& oldIs() const
        {
            return oldIs_;
        }

        //- Return number of stored fields
        inline label nOld() const
        {
            return nOld_;
        }

        //- Return stored delta indexes
        inline const labelList& deltaIs() const
        {
            return deltaIs_;
        }

        //- Return number of stored deltas
        inline label nDelta() const
        {
            return nDelta_;
        }

        //- Return current step
        inline label step() const
        {
            return stepi_;
        }

        //- Return local step fraction
        scalar f() const;

        //- Return start fraction of the current step
        scalar f0() const;

        //- Is the current step the last
        bool finalStep() const;

        //- Return current step VbyV0
        const volScalarField::Internal& V0byV() const
        {
            return V0byV_;
        }

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
        oldFields                                               \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        );                                                      \
                                                                \
        inline const PtrList<FieldName(Geo, Type)>&             \
        deltaFields                                             \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        ) const;                                                \
                                                                \
        inline PtrList<FieldName(Geo, Type)>&                   \
        deltaFields                                             \
        (                                                       \
            const FieldName(Geo, Type)&                         \
        );
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

#include "timeIntegratorI.H"

#ifdef NoRepository
#include "timeIntegratorTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
