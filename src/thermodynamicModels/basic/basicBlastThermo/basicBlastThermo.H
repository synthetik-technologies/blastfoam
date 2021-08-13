/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Added simple class to derive fluid and solid thermos
                            from
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
    Foam::basicBlastThermo

Description
    General thermodynamic class for a fluid/solid

SourceFiles
    basicBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicBlastThermo_H
#define basicBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "wordIOList.H"
#include "integrationSystem.H"
#include "regIOobject.H"
#include "speciesTable.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class basicBlastThermo
:
    public integrationSystem
{
protected:

    // Protected data

        // Thermo dictionary
        dictionary thermoDict_;

        //- Name of master thermo class
        word masterName_;

        //- Names of phase
        word name_;

        //- Phase density
        volScalarField& rho_;

        //- Reference to temperature
        volScalarField& T_;

        //- Reference to internal energy
        volScalarField& e_;

        //- Residual volume fraction
        dimensionedScalar residualAlpha_;

        //- Residual density
        dimensionedScalar residualRho_;

        //- Generic lookup for thermodynamics package thermoTypeName
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoTypeDict,
            Table* tablePtr,
            const int nCmpts,
            const char* cmptNames[],
            const word& thermoTypeName
        );

        //- Generic lookup for each of the related thermodynamics packages
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoDict,
            Table* tablePtr
        );

        //- Generic lookup for each of the related thermodynamics packages
        template<class Thermo, class Table>
        static typename Table::iterator lookupCstrIter
        (
            const dictionary& thermoDict,
            Table* tablePtr,
            const word& state
        );

        //- Return a subList
        UIndirectList<scalar> cellSetScalarList
        (
            const volScalarField& psi,
            const labelList& cells
        ) const;

public:

    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            basicBlastThermo,
            dictionary,
            (
                const word& phaseName,
                volScalarField& rho,
                volScalarField& e,
                volScalarField& T,
                const dictionary& dict,
                const word& masterName
            ),
            (phaseName, rho, e, T, dict, masterName)
        );


    //- Runtime type information
    TypeName("basicBlastThermo");


    // Constructor

        //- Construct with input field
        basicBlastThermo
        (
            const word& phaseName,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const word& masterName = word::null
        );

        //- Initialize models
        virtual void initializeModels()
        {}


    //- Destructor
    virtual ~basicBlastThermo();

    // Selectors

        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const word& phaseName,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const word& masterName = word::null
        );

        static autoPtr<basicBlastThermo> New
        (
            const word& phaseName,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const word& masterName = word::null
        );


    // Static public functions

        word phasePropertyName(const word& n) const
        {
            return IOobject::groupName(n, name_);
        }

        //- Return the thermo type name given a dictionary
        static word readThermoType(const dictionary& dict);

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName
        );

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName,
            const label nCmpts
        );

        //- Lookup or construct field
        static volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const word& name,
            const IOobject::readOption rOpt,
            const IOobject::writeOption wOpt,
            const dimensionSet& dims,
            const bool allowNoGroup = false
        );

        virtual word thermoName() const = 0;

    // Access functions
        const word& name() const
        {
            return name_;
        }

        const dimensionedScalar& residualAlpha() const
        {
            return residualAlpha_;
        }

        const dimensionedScalar& residualRho() const
        {
            return residualRho_;
        }

        //- Return constant reference to density field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return non-constant reference to density field
        volScalarField& rho()
        {
            return rho_;
        }

        //- Return constant reference to density field
        const volScalarField& he() const
        {
            return e_;
        }

        //- Return non-constant reference to density field
        volScalarField& he()
        {
            return e_;
        }

        //- Return constant reference to density field
        const volScalarField& T() const
        {
            return T_;
        }

        //- Return non-constant reference to density field
        volScalarField& T()
        {
            return T_;
        }

        //- Return mass fraction of a specie
        virtual const volScalarField& Y(const word& name) const
        {
            NotImplemented;
            return T_;
        }

        //- Access mass fraction of a specie
        virtual volScalarField& Y(const word& name)
        {
            NotImplemented;
            return T_;
        }

        //- Return mass fraction of a specie
        virtual const volScalarField& Y(const label i) const
        {
            NotImplemented;
            return rho_;
        }

        //- Access mass fraction of a specie
        virtual volScalarField& Y(const label i)
        {
            NotImplemented;
            return T_;
        }


    // Member Functions

        //- Correct thermodynamic fields
        virtual void correct() = 0;

        //- Names of species
        virtual speciesTable species() const
        {
            return speciesTable();
        }

        //- Does the mixture include this specie?
        virtual bool contains(const word& specieName) const
        {
            return false;
        }

        //- Return energy source
        virtual tmp<volScalarField> ESource() const = 0;

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, const volScalarField& delta);


        //- Thermodynamic and transport functions

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> he
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> he
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Sensible enthalpy [J/kg/K]
            virtual tmp<volScalarField> hs() const = 0;

            //- Sensible enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> hs
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Sensible enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> hs
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Absolute enthalpy [J/kg/K]
            virtual tmp<volScalarField> ha() const = 0;

            //- Absolute enthalpy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> ha
            (
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy for patch [J/kg/K]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Absolute enthalpy for cell-set [J/kg]
            virtual tmp<scalarField> ha
            (
                const scalarField& T,
                const labelList& cells
            ) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual tmp<volScalarField> hc() const = 0;

            //- Flame temperature [K]
            virtual tmp<volScalarField> flameT() const = 0;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE() const = 0;

            //- Temperature from enthalpy/internal energy
            virtual tmp<volScalarField> THE
            (
                const volScalarField& h,
                const volScalarField& p,
                const volScalarField& T0    // starting temperature
            ) const = 0;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const labelList& cells
            ) const = 0;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE
            (
                const scalarField& he,
                const scalarField& T0,      // starting temperature
                const label patchi
            ) const = 0;

            //- Temperature from enthalpy/internal energy for celli
            virtual scalar THEi
            (
                const scalar he,
                const scalar T0,
                const label celli
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp() const = 0;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant pressure for a cell [J/kg/K]
            virtual scalar Cpi(const label celli) const = 0;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const = 0;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

            //- Heat capacity at constant volume for a cell [J/kg/K]
            virtual scalar Cvi(const label celli) const = 0;

            //- Heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> Cpv() const = 0;

            //- Heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> Cpv
            (
                const scalarField& T,
                const label patchi
            ) const = 0;

             //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const = 0;

            //- Molecular weight for patch [kg/kmol]
            virtual tmp<scalarField> W(const label patchi) const = 0;

            //- Molecular weight for cell [kg/kmol]
            virtual scalar Wi(const label celli) const = 0;


    // Transport variables

        //- Return the thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volScalarField> alpha() const = 0;

        //- Read residual values
        virtual void read(const dictionary& dict)
        {
            this->residualRho_.read(dict);
            this->residualAlpha_.read(dict);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicBlastThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
