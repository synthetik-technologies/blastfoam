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
    Foam::compressibleBlastSystem

Description
    Base class for a collection of equation of states using a shared pressure
    and velocity (5 equation model)

SourceFiles
    compressibleBlastSystem.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleBlastSystem_H
#define compressibleBlastSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "compressibleSystem.H"
#include "fluidBlastThermo.H"


namespace Foam
{
class blastRadiationModel;

/*---------------------------------------------------------------------------*\
                           Class compressibleBlastSystem Declaration
\*---------------------------------------------------------------------------*/

class compressibleBlastSystem
:
    public compressibleSystem
{
protected:
// Protected data

    //- Primitative variables

        //- ThermodynamicModel
        autoPtr<fluidBlastThermo> thermoPtr_;

        //- Total mass
        volScalarField& rho_;

        //- Pressure
        volScalarField& p_;

        //- Temperature
        volScalarField& T_;

        //- Internal energy
        volScalarField& e_;


    //- Fluxes

        //- Radiation model
        autoPtr<blastRadiationModel> radiation_;


    // Protected member function

        //- Allocate models
        virtual void setModels();

        //- Add source terms to momentum and energy equations
            virtual void addSources
            (
                volVectorField::Internal& rhoUSource,
                volScalarField::Internal& rhoESource
            ) const;


public:

    // Constructor
    compressibleBlastSystem
    (
        const dictionary& dict,
        const fvMesh& mesh,
        const word& thermoType
    );

    //- Destructor
    virtual ~compressibleBlastSystem();


    // Member functions

        //- Decode primative variables
        virtual void decode();

        //- Update fluxes
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();


    // Member Access Functions

        //- Const access to thermo model
        virtual const fluidBlastThermo& thermo() const
        {
            return thermoPtr_();
        }

        //- Non-const access to thermo model
        virtual fluidBlastThermo& thermo()
        {
            return thermoPtr_();
        }

        //- Const access to total mass
        virtual const volScalarField& rho() const
        {
            return rho_;
        }

        //- Const access to total mass
        virtual volScalarField& rho()
        {
            return rho_;
        }

        //- Return internal energy field
        virtual const volScalarField& he() const
        {
            return e_;
        }

        //- Access internal energy field
        virtual volScalarField& he()
        {
            return e_;
        }

        //- Return calculated temperature field
        virtual const volScalarField& T() const
        {
            return T_;
        }

        //- Const access to pressure
        virtual const volScalarField& p() const
        {
            return p_;
        }

        //- Return speed of sound
        virtual tmp<volScalarField> speedOfSound() const
        {
            return volScalarField::New
            (
                "speedOfSound",
                thermoPtr_->speedOfSound()
            );
        }

        //- Dummy write for regIOobject
        bool writeData(Ostream& os) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
