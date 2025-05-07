/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::surfaceTemperatureModel

Description
    Model for describing surface temperature of a phase

SourceFiles
    surfaceTemperatureModel.C
    surfaceTemperatureModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceTemperatureModel_H
#define surfaceTemperatureModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseModel.H"
#include "runTimeSelectionTables.H"
#include "phasePair.H"
#include "fluidBlastThermo.H"
#include "timeIntegrationSystem.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class surfaceTemperatureModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceTemperatureModel
:
    public timeIntegrationSystem
{
protected:

    //- Const reference to phase model
    const phaseModel& phase_;


public:

    //- Runtime type information
    TypeName("surfaceTemperatureModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            surfaceTemperatureModel,
            dictionary,
            (
                const dictionary& dict,
                const phaseModel& phase
            ),
            (dict, phase)
        );


    // Static data members


    // Constructors

        //- Construct froma dictionary and a phase model
        surfaceTemperatureModel
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~surfaceTemperatureModel();


    // Selectors

        static autoPtr<surfaceTemperatureModel> New
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    // Member Functions

        //- Return the phase model
        inline const phaseModel& phase() const
        {
            return phase_;
        }

        //- Return the surface temperature
        virtual tmp<volScalarField> Ts() const = 0;

        //- Update before each ode step
        virtual void update() = 0;

        //- Update after ode integration
        virtual void postUpdate() = 0;

        //- Solve sub-step stepi
        virtual void solve() = 0;
        
        //- Initialize Models
        virtual void initializeModels()
        {}
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
