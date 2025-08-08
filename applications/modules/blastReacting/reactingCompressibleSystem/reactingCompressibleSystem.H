/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2025
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
    Foam::reactingCompressibleSystem

Description
    Uses standard OpenFOAM thermodynamic classes to solve for a reacting
    system

SourceFiles
    reactingCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef reactingCompressibleSystem_H
#define reactingCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystem.H"
#include "psiMulticomponentThermo.H"
#include "fluidMulticomponentThermophysicalTransportModel.H"
#include "combustionModel.H"
#include "radiationModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class reactingCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class reactingCompressibleSystem
:
    public compressibleSystem
{
protected:
// Protected data

    //- Thermo
    autoPtr<psiMulticomponentThermo> thermo_;

    //- Primitative variables

        //- Total mass
        volScalarField rho_;

        //- Pressure
        volScalarField& p_;

        //- Temperature
        volScalarField& T_;

        //- Internal energy
        volScalarField& e_;

        //- Thermo turbulence model
        autoPtr<fluidMulticomponentThermophysicalTransportModel>
            thermophysicalTransport_;

        //- Chemistry
        autoPtr<combustionModel> reaction_;

        //- Radiation model
        autoPtr<radiationModel> radiation_;

        //- Solve combustion each ode step
        bool odeCombustion_;


    //- Allocate models
    virtual void setModels()
    {}


public:

    TypeName("reactingCompressibleSystem");

    // Constructor
    reactingCompressibleSystem
    (
        const dictionary& dict,
        const fvMesh& mesh
    );

    //- Destructor
    virtual ~reactingCompressibleSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();


    // Member Access Functions

        //- Const access to thermo model
        const fluidThermo& thermo() const
        {
            return thermo_();
        }

        //- Non-const access to thermo model
        fluidThermo& thermo()
        {
            return thermo_();
        }

        //- Const access to total mass
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Access to total mass
        volScalarField& rho()
        {
            return rho_;
        }

        //- Const access to total mass
        const volScalarField& rhoEff() const
        {
            return rho_;
        }

        //- Access to total mass
        volScalarField& rhoEff()
        {
            return rho_;
        }

        //- Const access to internal energy
        const volScalarField& he() const
        {
            return e_;
        }

        //- Access to internal energy
        volScalarField& he()
        {
            return e_;
        }

        //- Return thermodynamic pressure
        const volScalarField& p() const
        {
            return p_;
        }

        //- Return temperature
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return speed of sound
        tmp<volScalarField> speedOfSound() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
