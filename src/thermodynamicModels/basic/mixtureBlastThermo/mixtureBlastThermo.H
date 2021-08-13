/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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
    Foam::mixtureBlastThermo

Description
    Templated class to handle multicomponent mixtures

SourceFiles
    mixtureBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureBlastThermo_H
#define mixtureBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicBlastThermo.H"
#include "speciesMixtureField.H"
#include "multicomponentBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixtureBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class ThermoType>
class mixtureBlastThermo
:
    public multicomponentBlastThermo<BasicThermo>
{
protected:
//- Protected data

    typedef ThermoType thermoType;

    //- Species data
    PtrList<ThermoType> speciesData_;

    //- Mixture
    speciesMixtureField<ThermoType> mixture_;



    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const PtrList<ThermoType>& constructSpeciesData(const dictionary&);

        //- Return a volScalarField of the given property
        template<class Method, class ... Args>
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a volScalarField of the given property for a specie
        template<class Method, class ... Args>
        tmp<volScalarField> volScalarFieldSpecieProperty
        (
            const word& s,
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a patch
        template<class Method, class ... Args>
        tmp<scalarField> patchFieldProperty
        (
            Method psiMethod,
            const label patchi,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a cell set
        template<class Method, class ... Args>
        tmp<scalarField> cellSetProperty
        (
            Method psiMethod,
            const labelList& cells,
            const Args& ... args
        ) const;


public:

    // Constructor
    mixtureBlastThermo
    (
        const word& phaseName,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const word& masterName = word::null
    );

    // Constructor
    mixtureBlastThermo
    (
        const HashPtrTable<ThermoType, word, string::hash>& thermoData,
        const word& phaseName,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const word& masterName = word::null
    );


    //- Destructor
    virtual ~mixtureBlastThermo();


    // Member functions

        //- Correct fields
        virtual void correct();

        //- Thermodynamic and transport functions

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy [J/kg/K]
            virtual tmp<volScalarField> hs() const;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Sensible enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Absolute enthalpy [J/kg/K]
            virtual tmp<volScalarField> ha() const;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const;

            //- Absolute enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Absolute enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const labelList& cells
            ) const;

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Flame temperature [K]
            virtual tmp<volScalarField> flameT() const;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE() const;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const;

            //- Temperature from enthalpy/internal energy for celli
            virtual scalar THEi
            (
                const scalar he,
                const scalar T0,
                const label celli
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for a cell [J/kg/K]
            virtual scalar Cpi(const label celli) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for a cell [J/kg/K]
            virtual scalar Cvi(const label celli) const;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const;

             //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const;

            //- Molecular weight for cell [kg/kmol]
            virtual scalar Wi(const label celli) const;


    // Fields derived from transport state variables

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> alpha() const;


        // Species properties

            //- Return molecular weight
            virtual tmp<volScalarField> specieW(const word& s) const;

            //- Return molecular weight
            virtual scalar specieWi(const word& s, const label celli) const;

            //- Mie Gruniesen coefficient
            virtual tmp<volScalarField> specieGamma(const word& s) const;

            //- Mie Gruniesen coefficient for a cell
            virtual scalar specieGammai(const word& s, const label celli) const;

            //- Heat capacity at constant pressure of mixture [J/kg/K]
            virtual tmp<volScalarField> specieCp(const word& s) const;

            //- Heat capacity at constant pressure of mixture for a cell [J/kg/K]
            virtual scalar specieCpi(const word& s, const label celli) const;

            //- Heat capacity at constant volume of mixture [J/kg/K]
            virtual tmp<volScalarField> specieCv(const word& s) const;

            //- Heat capacity at constant volume of mixture for a cell [J/kg/K]
            virtual scalar specieCvi(const word& s, const label celli) const;

            //- Specific heat ratio
            virtual tmp<volScalarField> specieCpByCv(const word& s) const;

            //- Heat of formation
            virtual tmp<volScalarField> specieHf(const word& s) const;

            //- Heat of formation
            virtual tmp<volScalarField> specieFlameT(const word& s) const;


    //- Return

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return ThermoType::typeName();
        }

        //- Return the table of thermos
        inline const PtrList<ThermoType>& speciesData() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtureBlastThermoI.H"

#ifdef NoRepository
    #include "mixtureBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
