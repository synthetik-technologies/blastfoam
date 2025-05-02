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

#include "multicomponentThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class basicSpecieBlastMixture Declaration
\*---------------------------------------------------------------------------*/

class basicSpecieBlastMixture
:
    public multicomponentThermo
{
protected:

    // Protected data

        //- Phase name
        const word phaseName_;

        //- Table of specie names
        speciesTable species_;

        //- The name of the default specie. The mass fraction of which is
        //  derived from the other species rather than solved. Also used as the
        //  carrier specie in multicomponent diffusion.
        word defaultSpecieName_;

        //- The index of the default specie
        label defaultSpeciei_;

        //- List of specie active flags
        List<bool> active_;

        //- Species mass fractions
        PtrList<volScalarField> Y_;


    // Protected Member Functions

        //- Normalise the mass fractions to sum to 1
        void normaliseMassFractions();

        //- Scale the mass fractions to sum to 1
        void correctMassFractions();

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

        // Species set

            //- The table of species
            virtual const speciesTable& species() const
            {
                return species_;
            }

            //- Does the mixture include this specie?
            inline bool containsSpecie(const word& specieName) const
            {
                return species_.found(specieName);
            }

            //- The index of the default specie
            virtual label defaultSpecie() const
            {
                return defaultSpeciei_;
            }


        // Active/Inactive handling

            //- Access the specie active flags
            virtual const List<bool>& speciesActive() const
            {
                return active_;
            }


        // Mass fractions

            //- Access the mass-fraction fields
            virtual PtrList<volScalarField>& Y()
            {
                return Y_;
            }

            //- Access the mass-fraction fields
            virtual const PtrList<volScalarField>& Y() const
            {
                return Y_;
            }

        // Specie thermo properties

            //- Return the pressure for cell [Pa]
            virtual scalar cellp(const label celli) const = 0;

            //- Pressure [kg/m/s^2]
            virtual scalar pi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> pi
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Return derivative of pressure w.r.t. density
            virtual scalar dpdRhoi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Return derivative of pressure w.r.t. Temperature
            virtual scalar dpdTi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            using multicomponentThermo::Cpi;
            virtual scalar Cpi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cpi
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            using multicomponentThermo::hei;
            virtual scalar hei
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> hei
            (
                const label speciei,
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> hei
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            using multicomponentThermo::hsi;
            virtual scalar hsi
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> hsi
            (
                const label speciei,
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> hsi
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            using multicomponentThermo::hai;
            virtual scalar hai
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> hai
            (
                const label speciei,
                const scalarField& rho,
                const scalarField& e,
                const scalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> hai
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;


        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mui
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mui
            (
                const label speciei,
                const volScalarField& rho,
                const volScalarField& e,
                const volScalarField& T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            using multicomponentThermo::kappai;
            virtual scalar kappai
            (
                const label speciei,
                const scalar rho,
                const scalar e,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappai
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
