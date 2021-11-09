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

    //- FvModels
    mutable UautoPtr<fvModels> modelsPtr_;

    //- fvConstraints
    mutable UautoPtr<const fvConstraints> constraintsPtr_;

    //- List of fields that have models or constraints
    mutable hashedWordList solveFields_;

    //- Stored old fields
    mutable HashPtrTable<PtrList<volScalarField>> oldScalarFields_;
    mutable HashPtrTable<PtrList<volVectorField>> oldVectorFields_;
    mutable HashPtrTable<PtrList<volSphericalTensorField>>
        oldSphTensorFields_;
    mutable HashPtrTable<PtrList<volSymmTensorField>> oldSymmTensorFields_;
    mutable HashPtrTable<PtrList<volTensorField>> oldTensorFields_;

    //- Stored delta fields
    mutable HashPtrTable<PtrList<volScalarField>> deltaScalarFields_;
    mutable HashPtrTable<PtrList<volVectorField>> deltaVectorFields_;
    mutable HashPtrTable<PtrList<volSphericalTensorField>>
        deltaSphTensorFields_;
    mutable HashPtrTable<PtrList<volSymmTensorField>> deltaSymmTensorFields_;
    mutable HashPtrTable<PtrList<volTensorField>> deltaTensorFields_;

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
            void clearFields(HashPtrTable<PtrList<FieldType>>& table);


public:

    //- Runtime type information
    TypeName("timeIntegrator");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            timeIntegrator,
            dictionary,
            (const fvMesh& mesh, const label nSteps),
            (mesh, nSteps)
        );

    // Constructor
    timeIntegrator(const fvMesh& mesh, const label nSteps);


    //- Destructor
    virtual ~timeIntegrator();


    // Selector
    static autoPtr<timeIntegrator> New(const fvMesh& mesh);


    // Member Functions

        //- Add a system to the
        void addSystem(timeIntegrationSystem& system);

        //- Update anything that needs to be done before mesh updates
        void preUpdateMesh();

        //- Integrate fluxes in time
        virtual void integrate();

        //- Clear all ODE fields (old/delta)
        void clear();

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

        //- Add old field to saved fields
        inline void addOldField(const volScalarField& f) const;
        inline void addOldField(const volVectorField& f) const;
        inline void addOldField(const volSphericalTensorField& f) const;
        inline void addOldField(const volSymmTensorField& f) const;
        inline void addOldField(const volTensorField& f) const;

        //- Add old field to saved fields
        inline void addDeltaField(const volScalarField& f) const;
        inline void addDeltaField(const volVectorField& f) const;
        inline void addDeltaField(const volSphericalTensorField& f) const;
        inline void addDeltaField(const volSymmTensorField& f) const;
        inline void addDeltaField(const volTensorField& f) const;

        //- Lookup old fields (const)
        inline const PtrList<volScalarField>& oldFields
        (
            const volScalarField& f
        ) const;
        inline const PtrList<volVectorField>& oldFields
        (
            const volVectorField& f
        ) const;
        inline const PtrList<volSphericalTensorField>& oldFields
        (
            const volSphericalTensorField& f
        ) const;
        inline const PtrList<volSymmTensorField>& oldFields
        (
            const volSymmTensorField& f
        ) const;
        inline const PtrList<volTensorField>& oldFields
        (
            const volTensorField& f
        ) const;

        //- Lookup old fields (non-const)
        inline PtrList<volScalarField>& oldFields(const volScalarField& f);
        inline PtrList<volVectorField>& oldFields(const volVectorField& f);
        inline PtrList<volSphericalTensorField>& oldFields
        (
            const volSphericalTensorField& f
        );
        inline PtrList<volSymmTensorField>& oldFields
        (
            const volSymmTensorField& f
        );
        inline PtrList<volTensorField>& oldFields(const volTensorField& f);

        //- Lookup delta fields (const)
        inline const PtrList<volScalarField>& deltaFields
        (
            const volScalarField& f
        ) const;
        inline const PtrList<volVectorField>& deltaFields
        (
            const volVectorField& f
        ) const;
        inline const PtrList<volSphericalTensorField>& deltaFields
        (
            const volSphericalTensorField& f
        ) const;
        inline const PtrList<volSymmTensorField>& deltaFields
        (
            const volSymmTensorField& f
        ) const;
        inline const PtrList<volTensorField>& deltaFields
        (
            const volTensorField& f
        ) const;

        //- Lookup delta fields (non-const)
        inline PtrList<volScalarField>& deltaFields(const volScalarField& f);
        inline PtrList<volVectorField>& deltaFields(const volVectorField& f);
        inline PtrList<volSphericalTensorField>& deltaFields
        (
            const volSphericalTensorField& f
        );
        inline PtrList<volSymmTensorField>& deltaFields
        (
            const volSymmTensorField& f
        );
        inline PtrList<volTensorField>& deltaFields(const volTensorField& f);


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
