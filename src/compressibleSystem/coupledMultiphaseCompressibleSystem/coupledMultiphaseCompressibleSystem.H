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
    Foam::coupledMultiphaseCompressibleSystem

Description
    Compressible system containing more than two equation of states.

SourceFiles
    coupledMultiphaseCompressibleSystem.C

\*---------------------------------------------------------------------------*/

#ifndef coupledMultiphaseCompressibleSystem_H
#define coupledMultiphaseCompressibleSystem_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multiphaseCompressibleSystem.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class coupledMultiphaseCompressibleSystem Declaration
\*---------------------------------------------------------------------------*/

class coupledMultiphaseCompressibleSystem
:
    public multiphaseCompressibleSystem
{
    // Private data

        //- Total volume fraction
        volScalarField volumeFraction_;

        //- Pointer to dispersed phase volume fraction
        const volScalarField* alphadPtr_;

        //- Phasic density
        volScalarField alphaRho_;

        //- Calculate new alpha and rho fields
        virtual void calcAlphaAndRho();


public:

    TypeName("coupledMultiphaseCompressibleSystem");

    // Constructor
    coupledMultiphaseCompressibleSystem(const fvMesh& mesh, const dictionary& dict);


    //- Destructor
    ~coupledMultiphaseCompressibleSystem();


    // Member functions

        //- Set a pointer to the dispersed phase volume fraction
        void setDispersedVolumeFraction(const volScalarField& alphad)
        {
            alphadPtr_ = &alphad;
            decode();
        }

        //- Return const access to volume fraction
        const volScalarField& volumeFraction() const
        {
            return volumeFraction_;
        }

        //- Non-cont access to volume fraction
        volScalarField& volumeFraction()
        {
            return volumeFraction_;
        }

        //- Decode primative variables
        virtual void decode();

        //- Encode conserved variables
        virtual void encode();

        //- Solve a sub-step
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );


    //- Return functions for transport model

        //- Heat capacit at contant volume of mixture [J/kg/K]
        virtual tmp<volScalarField> Cv() const;

        //- Dynamic viscosity of mixture [kg/m/s]
        virtual tmp<volScalarField> mu() const;

        //- Dynamic viscosity of mixture for patch [kg/m/s]
        virtual tmp<scalarField> mu(const label patchi) const;

        //- Kinematic viscosity of mixture [m^2/s]
        virtual tmp<volScalarField> nu() const;

        //- Kinematic viscosity of mixture for patch [m^2/s]
        virtual tmp<scalarField> nu(const label patchi) const;

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volScalarField> alpha() const;

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff(const volScalarField&) const;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Thermal diffusivity for energy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const;

        //- Thermal diffusivity for temperature of mixture [W/m/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappa(const label patchi) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff(const volScalarField&) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
