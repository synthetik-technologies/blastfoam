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

#include "timeIntegrationSystemBase.H"
#include "fvMesh.H"
#include "Time.H"
#include "fvTimeIntegrator.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class timeIntegrationSystem Declaration
\*---------------------------------------------------------------------------*/

class timeIntegrationSystem
:
    public timeIntegrationSystemBase
{
protected:
// Protected data

        //- Constant reference to mesh
        UautoPtr<const fvMesh> meshPtr_;

        //- Stored reference to time integrator
        UautoPtr<const fvTimeIntegrator> fvTimeInt_;

        //- Store old times when starting integration loop
        bool storeOld_;

        // Storage for fields

            using timeIntegrationSystemBase::storeOld;
            using timeIntegrationSystemBase::storeDelta;

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


        // General functions

            using timeIntegrationSystemBase::blendOld;
            using timeIntegrationSystemBase::blendDelta;
            using timeIntegrationSystemBase::storeAndBlendOld;
            using timeIntegrationSystemBase::storeAndBlendDelta;
            using timeIntegrationSystemBase::calcDelta;
            using timeIntegrationSystemBase::calcAndStoreDelta;

            //- Combine old fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void blendOld
            (
                FieldType& f,
                const bool conservative = true
            ) const;

            //- Combine delta fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void blendDelta(FieldType& f) const;

            //- Combine old fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void storeAndBlendOld
            (
                FieldType& f,
                const bool conservative = true
            );

            //- Combine delta fields using time step coefficients, looking up
            //  the saved fields
            template<class FieldType>
            void storeAndBlendDelta(FieldType& f);

            //- Calculate the sub delta given the actual change, looking up
            //  the saved fields
            template<class FieldType>
            tmp<FieldType> calcDelta(const FieldType& f) const;

            //- Calculate the sub delta given the actual change, looking up
            //  the saved fields
            template<class FieldType>
            tmp<FieldType> calcAndStoreDelta(const FieldType& f);

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
    timeIntegrationSystem(const word& name);
    timeIntegrationSystem(const word& name, const fvMesh& mesh);

    //- Set pointer to time integration system
    using timeIntegrationSystemBase::set;
    virtual void set(const fvMesh& mesh);
    virtual void set(const timeIntegrator& integrator);


    //- Destructor
    virtual ~timeIntegrationSystem();


    // Member functions

        //- Return the mesh
        const fvMesh& mesh() const
        {
            return meshPtr_();
        }

        //- Optional update before changes to the mesh
        virtual void preUpdateMesh()
        {}

        //- Return the fvModels
        inline const fvModels& models() const
        {
            return fvTimeInt_->models();
        }

        //- Access the fvModels
        inline fvModels& modelsRef() const
        {
            return fvTimeInt_->modelsRef();
        }

        //- Return the fvConstraints
        inline const fvConstraints& constraints() const
        {
            return fvTimeInt_->constraints();
        }

        //- Is a field modified by fvModels or fvConstraints
        inline bool needSolve(const word& field) const
        {
            return fvTimeInt_->needSolve(field);
        }
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
