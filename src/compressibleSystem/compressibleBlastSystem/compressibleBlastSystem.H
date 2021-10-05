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


        //- External energy source matrix
        tmp<fvScalarMatrix> extESource_;

        //- Drag matrix
        tmp<fvVectorMatrix> dragSource_;

        //- Allocate models
        virtual void setModels();


public:

    // Constructor
    compressibleBlastSystem(const label, const fvMesh& mesh);

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

        //- Return energy source term
        virtual tmp<volScalarField> ESource() const
        {
            return thermoPtr_->ESource();
        }

        //- Add external energy source
        tmp<fvScalarMatrix>& eSource()
        {
            return extESource_;
        }

        //- Add external energy source
        void addECoeff(const volScalarField::Internal& extEsrc);

        //- Add external energy source
        void addESource(const volScalarField::Internal& extEsrc);

        //- Add external energy source
        tmp<fvVectorMatrix>& dragSource()
        {
            return dragSource_;
        }

        //- Add to U coefficient
        void addUCoeff(const volScalarField::Internal& coeff);

        //- Add to U source
        void addUSource(const volVectorField::Internal& src);


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
