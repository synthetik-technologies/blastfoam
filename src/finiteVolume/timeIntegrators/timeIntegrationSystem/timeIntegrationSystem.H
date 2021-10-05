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
    Foam::timeIntegrationSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    timeIntegrationSystem.C

\*---------------------------------------------------------------------------*/

#ifndef timeIntegrationSystem_H
#define timeIntegrationSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UautoPtr.H"
#include "fvMesh.H"
#include "Time.H"
#include "timeIntegrator.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class timeIntegrationSystem Declaration
\*---------------------------------------------------------------------------*/

class timeIntegrationSystem
{
protected:
// Protected data

        //- Constant reference to mesh
        UautoPtr<const fvMesh> meshPtr_;

        //- Name of system
        word name_;

        //- Stored reference to time integrator
        UautoPtr<timeIntegrator> timeInt_;

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

        // Storage for fields

            //- Store old fields
            template<class FieldType>
            void storeOld
            (
                FieldType& f,
                PtrList<FieldType>& fList,
                const bool conservative = true
            );

            //- Store delta fields
            template<class FieldType>
            void storeDelta
            (
                const FieldType& f,
                PtrList<FieldType>& fList
            );

            //- Store old fields
            template<class FieldType>
            void storeOld
            (
                FieldType& f,
                const bool conservative = true
            );

            //- Store delta fields
            template<class FieldType>
            void storeDelta(const FieldType& f);


        // storage for single quantities

            //- Store old fields
            //  moving keyword used for template function
            template<class Type>
            void storeOld(Type& f, List<Type>& fList, bool conservative = true);

            //- Store delta fields
            template<class Type>
            void storeDelta(const Type& f, List<Type>& fList);


        // General functions

            //- Combine old fields using time step coefficients
            template<template<class> class ListType, class Type>
            void blendOld
            (
                Type& f,
                const ListType<Type>& fList
            ) const;

            //- Combine old fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void blendOld(FieldType& f) const;

            //- Combine delta fields using time step coefficients
            template<template<class> class ListType, class Type>
            void blendDelta
            (
                Type& f,
                const ListType<Type>& fList
            ) const;

            //- Combine delta fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void blendDelta(FieldType& f) const;

            //- Combine old fields using time step coefficients
            template<template<class> class ListType, class Type>
            void storeAndBlendOld
            (
                Type& f,
                ListType<Type>& fList,
                const bool conservative = true
            );

            //- Combine old fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void storeAndBlendOld(FieldType& f, const bool conservative = true);

            //- Combine delta fields using time step coefficients
            template<template<class> class ListType, class Type>
            void storeAndBlendDelta
            (
                Type& f,
                ListType<Type>& fList
            );

            //- Combine delta fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void storeAndBlendDelta(FieldType& f);

            //- Calculate the sub delta given the actual change
            template<template<class> class ListType, class Type>
            tmp<Type> calcDelta
            (
                const Type& f,
                const ListType<Type>& fList
            ) const;

            //- Calculate the sub delta given the actual change, looking up
            //  the saved fields
            template<class FieldType>
            tmp<FieldType> calcDelta(const FieldType& f) const;

            //- Calculate the sub delta given the actual change
            template<template<class> class ListType, class Type>
            tmp<Type> calcAndStoreDelta
            (
                const Type& f,
                ListType<Type>& fList
            );

            //- Calculate the sub delta given the actual change, looking up
            //  the saved fields
            template<class FieldType>
            tmp<FieldType> calcAndStoreDelta(const FieldType& f);

            //- Combine fields using time step coefficients
            //  f must be the most recent value of the field
            template<template<class> class ListType, class Type>
            void blendSteps
            (
                const labelList& indices,
                Type& f,
                const ListType<Type>& fList,
                const scalarList& scales
            ) const;

            //- Add old fields for a given variable
            template<class FieldType>
            void addOldField(const FieldType& f);

            //- Add delta fields for a given variable
            template<class FieldType>
            void addDeltaField(const FieldType& f);

            //- Add old and delta fields for a given variable
            template<class FieldType>
            void addOldDeltaField(const FieldType& f);


public:


    // Constructor
    timeIntegrationSystem();
    timeIntegrationSystem
    (
        const word& name,
        const fvMesh& mesh
    );
    void set(const word& name, const fvMesh& mesh);


    //- Destructor
    virtual ~timeIntegrationSystem();


    // Member functions

        //- Return name of system
        const word& name() const
        {
            return name_;
        }

        //- Return the mesh
        const fvMesh& mesh() const
        {
            return meshPtr_();
        }

        //- Update before each ode step
        virtual void update() = 0;

        //- Update after ode integration
        virtual void postUpdate() = 0;

        //- Solve sub-step stepi
        virtual void solve() = 0;

        //- Return the current step
        label step() const;

        //- Return old coefficients for the current step
        scalarList a() const;

        //- Return delta coefficients for the current step
        scalarList b() const;

        //- Return the time step fraction
        scalar f() const;

        //- Return the initial time fraction
        scalar f0() const;

        //- Is the current step the last
        bool finalStep() const;

        //- Current time for substep
        dimensionedScalar time() const;

        //- Current time step size
        dimensionedScalar deltaT() const;

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeIntegrationSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
