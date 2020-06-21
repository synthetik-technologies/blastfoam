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
    Foam::basicThermoModel

Description
    General class for a fluid/solid

SourceFiles
    basicThermoModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermoModel_H
#define basicThermoModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "wordIOList.H"
#include "regIOobject.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicThermoModel Declaration
\*---------------------------------------------------------------------------*/

class basicThermoModel
:
    public regIOobject
{
protected:

    // Protected data

        //- Is this the top level thermo model
        bool master_;

        //- Names of phase
        word name_;

        //- Pressure
        volScalarField& p_;

        //- Phase density
        volScalarField& rho_;

        //- Reference to temperature
        volScalarField& T_;

        //- Reference to internal energy
        volScalarField& e_;

        //- Thermal diffusivity [kg/m/s]
        volScalarField alpha_;

        //- Is the internal energy limited
        Switch limit_;


    //- Protected functions

        //- Correct e boundary values
        void eBoundaryCorrection();

        //- Lookup or construct field
        volScalarField& lookupOrConstruct
        (
            const fvMesh& mesh,
            const char* name,
            const IOobject::readOption rOpt,
            const IOobject::writeOption wOpt,
            const dimensionSet& dims
        ) const;

        //- Lookup or construct field
        volScalarField& lookupOrConstructE
        (
            const fvMesh& mesh,
            const char* name
        ) const;


public:

    //- Runtime type information
    TypeName("basicThermo");


    // Constructor

        //- Construct with input field
        basicThermoModel
        (
            const word& phaseName,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const bool master = true
        );

        //- Construct and add fields to db
        basicThermoModel
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict,
            const bool master = true
        );


    //- Destructor
    virtual ~basicThermoModel();


    // Static public functions

        //- Generic lookup for thermodynamics package thermoTypeName
        template<class Table>
        static typename Table::iterator lookupThermo
        (
            Table* tablePtr,
            const char* cmptNames[],
            const word& thermoTypeName
        );

        //- Generic lookup for each of the related thermodynamics packages
        template<class Table>
        static typename Table::iterator lookupThermo
        (
            const dictionary& thermoDict,
            Table* tablePtr
        );

        //- Return the thermo type name given a dictionary
        static word readThermoType(const dictionary& dict);

        //- Split name of thermo package into a list of the components names
        static wordList splitThermoName
        (
            const word& thermoName
        );

        //- Return the internal energy boundary types using an input
        //  temperature field
        static wordList eBoundaryTypes(const volScalarField& T);

        //- Return the internal energy boundary types using an input
        //  temperature field
        static wordList eBoundaryBaseTypes(const volScalarField& T);


    // Access functions
        const word& name() const
        {
            return name_;
        }

        //- Return constant reference to pressure field
        const volScalarField& p() const
        {
            return p_;
        }

        //- Return non-constant reference to pressure field
        volScalarField& p()
        {
            return p_;
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
        const volScalarField& e() const
        {
            return e_;
        }

        //- Return non-constant reference to density field
        volScalarField& e()
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


    // Member Functions

        //- Correct thermodynamic fields
        virtual void correct();

        //- Is the internal energy limited
        bool limit() const
        {
            return limit_;
        }


    //- Thermodynamic and transport functions

        //- Calculate temperature
        virtual tmp<volScalarField> calcT() const = 0;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> TRhoE
        (
            const scalarField& T,
            const scalarField& e,
            const label patchi
        ) const = 0;

        //- Calculate internal energy for celli
        virtual scalar TRhoEi
        (
            const scalar& T,
            const scalar& e,
            const label celli
        ) const = 0;

        //- Calculate internal energy
        virtual tmp<volScalarField> E() const = 0;

        //- Calculate internal energy for patchi
        virtual tmp<volScalarField> e
        (
            const volScalarField& rho,
            const volScalarField& e,
            const volScalarField& T
        ) const = 0;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> e
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const = 0;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> e
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const labelList& faceCells
        ) const = 0;

        //- Return molecular weight
        virtual tmp<volScalarField> W() const = 0;

        //- Return molecular weight for patchi
        virtual tmp<scalarField> W(const label patchi) const = 0;

        //- Return molecular weight
        virtual scalar Wi(const label celli) const = 0;

        //- Mie Gruniesen coefficient
        virtual tmp<volScalarField> Gamma() const = 0;

        //- Mie Gruniesen coefficient for a patch
        virtual tmp<scalarField> Gamma(const label patchi) const = 0;

        //- Heat capacity at constant pressure of mixture [J/kg/K]
        virtual tmp<volScalarField> Cp() const = 0;

        //- Heat capacity at constant pressure of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cp(const label patchi) const = 0;

        //- Heat capacity at constant pressure of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cp
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const = 0;


        //- Heat capacity at constant pressure of mixture for a cell [J/kg/K]
        virtual scalar Cpi(const label celli) const = 0;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        virtual tmp<volScalarField> Cv() const = 0;

        //- Heat capacity at constant volume of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cv(const label patchi) const = 0;

        //- Heat capacity at constant volume of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cv
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const = 0;

        //- Heat capacity at constant volume of mixture for a cell [J/kg/K]
        virtual scalar Cvi(const label celli) const = 0;

        //- Specific heat ratio
        virtual tmp<volScalarField> CpByCv() const = 0;

        //- Specific heat ratio for patchi
        virtual tmp<scalarField> CpByCv(const label patchi) const = 0;

        //- Specific heat ratio for patchi
        virtual tmp<scalarField> CpByCv
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const = 0;

    // Transport variables

        //- Thermal diffusivity for temperature of mixture [W/m/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal diffusivity for temperature of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappa(const label patchi) const;

        //- Thermal diffusivity for temperature of mixture for cell [W/m/K]
        virtual scalar kappai(const label cell) const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alpha() const
        {
            return alpha_;
        }

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const
        {
            return alpha_.boundaryField()[patchi];
        }

        //- Thermal diffusivity for enthalpy of mixture [kg/m/s]
        virtual tmp<volScalarField> alphahe() const;

        //- Thermal diffusivity for enthalpy of mixture for patch [kg/m/s]
        virtual tmp<scalarField> alphahe(const label patchi) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff(const volScalarField& alphat) const;

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture for patch [W/m/K]
        virtual tmp<scalarField> kappaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff(const volScalarField&) const;

        //- Effective thermal turbulent diffusivity of mixture
        //  for patch [kg/m/s]
        virtual tmp<scalarField> alphaEff
        (
            const scalarField& alphat,
            const label patchi
        ) const;

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //- Dummy read for regIOobject
        virtual bool read()
        {
            return true;
        }

        //- Read residual values
        virtual void read(const dictionary& dict)
        {
            return;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicThermoModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
