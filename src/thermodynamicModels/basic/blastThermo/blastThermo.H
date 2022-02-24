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
    Foam::blastThermo

Description
    General fluid class that contains functions to return pressure and speed
    of sound

SourceFiles
    blastThermo.C
    newblastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef blastThermo_H
#define blastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrationSystem.H"
#include "basicThermo.H"
#include "OSspecific.H"

namespace Foam
{
class basicBlastThermo;

/*---------------------------------------------------------------------------*\
                           Class blastThermo Declaration
\*---------------------------------------------------------------------------*/

class blastThermo
:
    public timeIntegrationSystem,
    public basicThermo::implementation
{
protected:

    // Protected data

        //- Name of the phase
        const word phaseName_;

        //- Internal energy
        volScalarField e_;

        //- Density field [kg/m^3]
        volScalarField rho_;

        //- Heat capacity at constant pressure [J/kg/K]
        volScalarField Cp_;

        //- Heat capacity at constant volume [J/kg/K]
        volScalarField Cv_;

        //- Minimum temperature
        scalar TLow_;

        //- Residual volume fraction
        dimensionedScalar residualAlpha_;

        //- Residual density
        dimensionedScalar residualRho_;


        // Protected member functions

            //- Temperature [K]
            //  Non-const access allowed for transport equations
            volScalarField& TRef()
            {
                return this->T_;
            }

            //- Return enthalpy/Internal energy [J/kg]
            volScalarField& heRef()
            {
                return e_;
            }

            //- Heat capacity at constant pressure [J/kg/K]
            volScalarField& CpRef()
            {
                return Cp_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            volScalarField& CvRef()
            {
                return Cv_;
            }

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            volScalarField& alphaRef()
            {
                return alpha_;
            }


public:

    TypeName("blastThermo");

    // Public Classes
    blastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );

    //- Destructor
    virtual ~blastThermo();

    //- Initialize models
    virtual void initializeModels()
    {}

    // Selector

        template<class Thermo>
        static autoPtr<Thermo> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& phaseName,
            const word& masterName
        );


    // Static member functions

        // Directory to compile from
        static fileName templateDir()
        {
            return getEnv("BLAST_DIR")/"etc/codeTemplates";
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
        static UIndirectList<scalar> cellSetScalarList
        (
            const volScalarField& psi,
            const labelList& cells
        );

    // Access functions

        //- Energy Source
        virtual tmp<volScalarField> ESource() const = 0;

        //- Return residual values
        virtual dimensionedScalar residualAlpha() const
        {
            return residualAlpha_;
        }

        virtual dimensionedScalar residualRho() const
        {
            return residualRho_;
        }

        virtual dimensionedScalar residualAlphaRho() const
        {
            return residualAlpha()*residualRho();
        }

        // Access functions

            const word& phaseName() const
            {
                return phaseName_;
            }

            const word& masterName() const
            {
                return phaseName();
            }

            //- Is the phase incompressible
            virtual bool incompressible() const
            {
                return false;
            }

            //- Constant volume
            virtual bool isochoric() const
            {
                return false;
            }


        // Fields derived from thermodynamic state variables

            using basicThermo::implementation::he;

            //- Return enthalpy/Internal energy [J/kg]
            virtual const volScalarField& he() const;

            //- Return enthalpy/Internal energy [J/kg]
            virtual volScalarField& he();

            //- Enthalpy/Internal energy for a cell [J/kg]
            virtual scalar cellHE
            (
                const scalar T,
                const label celli
            ) const = 0;

            //- Enthalpy/Internal energy for patch and face [J/kg]
            virtual scalar patchFaceHE
            (
                const scalar T,
                const label patchi,
                const label facei
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            using basicThermo::Cp;
            virtual tmp<volScalarField> Cp() const
            {
                return Cp_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            using basicThermo::Cv;
            virtual tmp<volScalarField> Cv() const
            {
                return Cv_;
            }

            //- Heat capacity at constant volume [J/kg/K]
            using basicThermo::Cpv;
            virtual tmp<volScalarField> Cpv() const
            {
                return Cv_;
            }

            using basicThermo::THE;

            //- Temperature from enthalpy/internal energy for celli
            virtual tmp<volScalarField> THE() const = 0;

            //- Temperature from enthalpy/internal energy for celli
            virtual scalar cellTHE
            (
                const scalar he,
                const scalar T0,
                const label celli
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Density for patch [kg/m^3]
            virtual tmp<scalarField> rho(const label patchi) const;

            //- Density for cell [kg/m^3]
            virtual scalar cellrho(const label celli) const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Old-time density [kg/m^3]
            virtual tmp<volScalarField> rho0() const;

            //- Recompute the density
            virtual void updateRho() = 0;

            //- Gamma = Cp/Cv []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = Cp/Cv for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Flame temperature [k]
            virtual tmp<volScalarField> flameT() const = 0;

            //- Return molecular weight
            virtual scalar cellW(const label) const = 0;

            //- Heat capacity at constant pressure of mixture for a cell [J/kg/K]
            virtual scalar cellCp(const scalar, const label) const = 0;

            //- Heat capacity at constant volume of mixture for a cell [J/kg/K]
            virtual scalar cellCv(const scalar, const label) const = 0;

            //- Heat capacity at constant pressure/volume for a cell [J/kg/K]
            virtual scalar cellCpv(const scalar, const label) const = 0;

            //- Heat capacity at constant pressure/volume for a patch and face [J/kg/K]
            virtual scalar patchFaceCpv
            (
                const scalar,
                const label,
                const label
            ) const = 0;


        // Fields derived from transport state variables

            //- Kinematic viscosity [m^2/s]
            virtual tmp<volScalarField> nu() const = 0;

            //- Thermal diffusivity for temperature of mixture [W/m/K]
            virtual tmp<volScalarField> kappa() const;

            //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappa(const label patchi) const;

            //- Thermal diffusivity for temperature of mixture for cell [W/m/K]
            virtual scalar cellkappa(const label cell) const;

            //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphahe() const;

            //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphahe(const label patchi) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture [W/m/K]
            virtual tmp<volScalarField> kappaEff
            (
                const volScalarField&
            ) const;

            //- Effective thermal turbulent diffusivity for temperature
            //  of mixture for patch [W/m/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal turbulent diffusivity of mixture
            //  for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

        //- Read residual values
        virtual bool read();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blastThermoTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
