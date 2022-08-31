/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2021 OpenFOAM Foundation
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
    Foam::basicSpecieBlastMixture

Description
    Extension of basicSpeciesMixture for density based thermos

SourceFiles
    basicSpecieBlastMixture.C

\*---------------------------------------------------------------------------*/

#ifndef basicSpecieBlastMixture_H
#define basicSpecieBlastMixture_H

#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class basicSpecieBlastMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSpecieBlastMixture
:
    public basicSpecieMixture
{
public:

    //- Run time type information
    TypeName("basicSpecieBlastMixture");

    //- The base class of the mixture
    typedef basicSpecieBlastMixture basicMixtureType;


    // Constructors

        //- Construct from dictionary, species names, mesh and phase name
        basicSpecieBlastMixture
        (
            const dictionary&,
            const wordList& specieNames,
            const fvMesh&,
            const word&
        );


    //- Destructor
    virtual ~basicSpecieBlastMixture()
    {}


    // Member Functions

        // Specie thermo properties

            //- Pressure [kg/m/s^2]
            virtual scalar p
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> p
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Return derivative of pressure w.r.t. density
            virtual scalar dpdRho
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return derivative of pressure w.r.t. Temperature
            virtual scalar dpdT
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            using basicSpecieMixture::Cp;
            virtual scalar Cp
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            using basicSpecieMixture::HE;
            virtual scalar HE
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> HE
            (
                const label speciei,
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> HE
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            using basicSpecieMixture::Hs;
            virtual scalar Hs
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> Hs
            (
                const label speciei,
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> Hs
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            using basicSpecieMixture::Ha;
            virtual scalar Ha
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> Ha
            (
                const label speciei,
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> Ha
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;


        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            using basicSpecieMixture::mu;
            virtual scalar mu
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mu
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            using basicSpecieMixture::kappa;
            virtual scalar kappa
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
