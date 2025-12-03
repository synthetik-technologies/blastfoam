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
    Foam::timeIntegrationSystemBase

Description
    Base class for a soling time dependent ODEs

SourceFiles
    timeIntegrationSystemBase.C

\*---------------------------------------------------------------------------*/

#ifndef timeIntegrationSystemBase_H
#define timeIntegrationSystemBase_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "UautoPtr.H"
#include "timeIntegrator.H"
#include "dimensionedTypes.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class timeIntegrationSystemBase Declaration
\*---------------------------------------------------------------------------*/

class timeIntegrationSystemBase
{
protected:
// Protected data

        //- Name of system
        word name_;

        //- Stored reference to time integrator
        UautoPtr<const timeIntegrator> timeInt_;


    // Protected member functions

        // Storage for single quantities

            //- Store old fields
            //  moving keyword used for template function
            template<class Type>
            void storeOld
            (
                const Type& f,
                List<Type>& fList,
                bool conservative = true
            );

            //- Store delta fields
            template<class Type>
            void storeDelta(const Type& f, List<Type>& fList);

        // Storage for fields

            //- Store old fields
            template<class FieldType>
            void storeOld
            (
                const FieldType& f,
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

        // General functions for blending

            //- Combine old fields using time step coefficients
            template<class Type>
            void blendOld
            (
                Type& f,
                const List<Type>& fList,
                const bool conservative = true
            ) const;

            //- Combine old fields using time step coefficients
            template<class FieldType>
            void blendOld
            (
                FieldType& f,
                const PtrList<FieldType>& fList,
                const bool conservative = true
            ) const;

            //- Combine delta fields using time step coefficients
            template<template<class> class ListType, class Type>
            void blendDelta
            (
                Type& f,
                const ListType<Type>& fList
            ) const;

            //- Combine old fields using time step coefficients
            template<template<class> class ListType, class Type>
            void storeAndBlendOld
            (
                Type& f,
                ListType<Type>& fList,
                const bool conservative = true
            );

            //- Combine delta fields using time step coefficients
            template<template<class> class ListType, class Type>
            void storeAndBlendDelta
            (
                Type& f,
                ListType<Type>& fList
            );

            //- Calculate the sub delta given the actual change
            template<template<class> class ListType, class Type>
            tmp<Type> calcDelta
            (
                const Type& f,
                const ListType<Type>& fList
            ) const;

            //- Calculate the sub delta given the actual change
            template<template<class> class ListType, class Type>
            void calcDelta
            (
                const Type& f,
                const ListType<Type>& fList,
                Type& fLimited
            ) const;

            //- Calculate the sub delta given the actual change
            template<template<class> class ListType, class Type>
            tmp<Type> calcAndStoreDelta
            (
                const Type& f,
                ListType<Type>& fList
            );

            //- Calculate the sub delta given the actual change
            template<template<class> class ListType, class Type>
            void calcAndStoreDelta
            (
                const Type& f,
                ListType<Type>& fList,
                Type& fLimited
            );


public:


    // Constructors

        // Construct from name and object registry
        timeIntegrationSystemBase(const word& name);

        // Construct from name and object registry
        timeIntegrationSystemBase
        (
            const word& name,
            const objectRegistry& obr
        );

        // Construct from name and timeIntegrator
        timeIntegrationSystemBase
        (
            const word& name,
            const timeIntegrator& integrator
        );


    //- Destructor
    virtual ~timeIntegrationSystemBase();


    // Member functions

        //- Set pointer to time integration system
        virtual void set(const objectRegistry& obr);

        //- Set pointer to time integration system
        virtual void set(const timeIntegrator& integrator);

        //- Return name of system
        const word& name() const
        {
            return name_;
        }

        //- Update before mesh changes
        virtual void preUpdateMesh()
        {}

        //- Save before solving
        virtual void save()
        {}

        //- Update before ode integration
        virtual void preUpdate()
        {}

        //- Update before each ode step
        virtual void update() = 0;

        //- Solve sub-step stepi
        virtual void solve() = 0;

        //- Explicit update after ode integration
        virtual void postExplicit() = 0;

        //- Store explicit variables
        virtual void storeExplicit() = 0;

        //- Implicit update after ode integration
        virtual void postImplicit() = 0;

        //- Update after solving
        virtual void postUpdate()
        {}

        //- Optional clearing step
        virtual void clear()
        {}

        //- Optional resetting when restarting a time step
        virtual void reset()
        {}

        //- Return the current step
        label step() const;

        //- Return old coefficients for the current step
        const List<scalar>& a() const;

        //- Return delta coefficients for the current step
        const List<scalar>& b() const;

        //- Return the time step fraction
        scalar f() const;

        //- Return the initial time fraction
        scalar f0() const;

        //- Is the current step the last
        bool finalStep() const;

        //- Current time for substep
        dimensionedScalar t() const;

        //- Current time step size
        dimensionedScalar dt() const;

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "timeIntegrationSystemBaseTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
