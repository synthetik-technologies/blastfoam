/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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
    Foam::integrationSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    integrationSystem.C

\*---------------------------------------------------------------------------*/

#ifndef integrationSystem_H
#define integrationSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh.H"
#include "Time.H"


namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class integrationSystem Declaration
\*---------------------------------------------------------------------------*/

class integrationSystem
{
protected:
// Protected data

        //- Constant reference to mesh
        const fvMesh& mesh_;

        //- Name of system
        const word name_;

        //- Stored old field indexes
        labelList oldIs_;

        //- Number of stored fields
        label nOld_;

        //- Stored delta indexes
        labelList deltaIs_;

        //- Number of stored deltas
        label nDelta_;

        //- Store old fields
        template<class fieldType>
        void storeOld
        (
            const label stepi,
            fieldType& f,
            PtrList<fieldType>& fList,
            const bool moving = true
        ) const;

        //- Store delta fields
        template<class fieldType>
        void storeDelta
        (
            const label stepi,
            const fieldType& f,
            PtrList<fieldType>& fList
        ) const;

        //- Combine fields using time step coefficients
        //  f must be the most recent value of the field
        template<class fieldType>
        void blendSteps
        (
            const label stepi,
            const labelList& indices,
            fieldType& f,
            const PtrList<fieldType>& fList,
            const scalarList& scales
        ) const;

        //- Combine old fields using time step coefficients
        template<class fieldType>
        void blendOld
        (
            const label stepi,
            fieldType& f,
            const PtrList<fieldType>& fList,
            const scalarList& scales
        ) const;

        //- Combine delta fields using time step coefficients
        template<class fieldType>
        void blendDelta
        (
            const label stepi,
            fieldType& f,
            const PtrList<fieldType>& fList,
            const scalarList& scales
        ) const;

        //- Clear and resize stored old fields
        template<class fieldType>
        void clearOld(PtrList<fieldType>& fList) const;

        //- Clear and resize stored delta fields
        template<class fieldType>
        void clearDelta(PtrList<fieldType>& fList) const;

        //- Lookup global integrator and initialize
        void lookupAndInitialize(const word& name = "globalTimeIntegrator");


public:


    // Constructor
    integrationSystem
    (
        const word& name,
        const fvMesh& mesh
    );


    //- Destructor
    virtual ~integrationSystem();


    // Member functions

        //- Return name of system
        const word& name() const
        {
            return name_;
        }

        //- Set old lists and fluxes (initialization of fields)
        void setODEFields
        (
            const label nSteps,
            const boolList& storeFields,
            const boolList& storeDeltas
        );

        //- Update before each ode step
        virtual void update() = 0;

        //- Update after ode integration
        virtual void postUpdate() = 0;

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        ) = 0;

        //- Remove stored fields
        virtual void clearODEFields() = 0;

        //- Return the time step fraction
        scalar f(const label stepi, const scalarList& b) const;


        //- Dummy write for regIOobject
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "integrationSystemTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
