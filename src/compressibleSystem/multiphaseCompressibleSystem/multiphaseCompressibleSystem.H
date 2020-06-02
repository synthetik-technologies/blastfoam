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
    Foam::multiphaseCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    multiphaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseCompressibleSystem_H
#define multiphaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseCompressibleSystem.H"
#include "multiphaseFluidThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiphaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class multiphaseCompressibleSystem
:
    public phaseCompressibleSystem
{
protected:
// Protected data

    //- Primitative variables

        //- Equation of states
        multiphaseFluidThermo thermo_;

        //- Volume fractions
        PtrList<volScalarField>& alphas_;

        //- Phase densities
        UPtrList<volScalarField>& rhos_;


    //- Conserved variables

        //- Phase masses
        PtrList<volScalarField> alphaRhos_;


    //- Fluxes

        //- Volume fraction fluxes
        PtrList<surfaceScalarField> alphaPhis_;

        //- Phase mass fluxes
        PtrList<surfaceScalarField> alphaRhoPhis_;


    // ODE variables

        //- Old values for ode solver
        PtrList<PtrList<volScalarField>> alphasOld_;
        PtrList<PtrList<volScalarField>> alphaRhosOld_;

        //- Stored delta fields
        PtrList<PtrList<volScalarField>> deltaAlphas_;
        PtrList<PtrList<volScalarField>> deltaAlphaRhos_;

        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho();


public:

    TypeName("multiphaseCompressibleSystem");

    // Constructor
    multiphaseCompressibleSystem(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    ~multiphaseCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const;

        //- Solve a sub-step
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes
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
            return thermo_;
        }

        //- Non-const access to thermo model
        virtual fluidThermoModel& thermo()
        {
            return thermo_;
        }

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const
        {
            return thermo_.speedOfSound();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
