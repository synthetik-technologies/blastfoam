/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::phaseCompressibleSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    phaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseCompressibleSystem_H
#define singlePhaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseCompressibleSystem.H"
#include "fluidThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class singlePhaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseCompressibleSystem
:
    public phaseCompressibleSystem
{
protected:
// Protected data

    //- Equation of state
    autoPtr<fluidThermoModel> thermo_;



    // ODE variables

        //- Old values for ode solver
        PtrList<volScalarField> rhoOld_;

        //- Stored delta fields
        PtrList<volScalarField> deltaRho_;

        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho()
        {
            return;
        }

public:

    TypeName("singlePhaseCompressibleSystem");

    // Constructor
    singlePhaseCompressibleSystem(const fvMesh& mesh, const dictionary& dict);

    //- Destructor
    virtual ~singlePhaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const;

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const boolList& storeFields,
            const boolList& storeDeltas
        );

        //- Remove stored fields
        virtual void clearODEFields();


    // Member Access Functions

        //- Const access to thermo model
        virtual const fluidThermoModel& thermo() const
        {
            return thermo_();
        }

        //- Non-const access to thermo model
        virtual fluidThermoModel& thermo()
        {
            return thermo_();
        }

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
