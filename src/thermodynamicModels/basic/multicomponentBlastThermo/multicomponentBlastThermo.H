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
    Foam::multicomponentBlastThermo

Description
    Class to handle species transport

SourceFiles
    multicomponentBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentBlastThermo_H
#define multicomponentBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "timeIntegrationSystem.H"
#include "PtrListDictionary.H"
#include "basicSpecieMixture.H"
#include "fvScalarMatrix.H"
#include "speciesTable.H"
#include "fvModels.H"
#include "fvConstraints.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multicomponentBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class multicomponentBlastThermo
:
    public basicSpecieMixture
{
//- Private data

    //- Name of the phase
    word phaseName_;

    //- Name of master thermo class
    word masterName_;


protected:
//- Protected data

    //- Mass transfer rates passed from an external source
    PtrListDictionary<volScalarField> massTransferRates_;

    //- Sources to solve implicitly
    PtrListDictionary<fvScalarMatrix> implicitSources_;

    //- Name of flux
    UautoPtr<const surfaceScalarField> alphaRhoPhiPtr_;

    //- Name of mass
    UautoPtr<const volScalarField> alphaRhoPtr_;

    class integrator
    :
        public timeIntegrationSystem
    {
    // Private data
        const fvMesh& mesh_;
        PtrList<volScalarField>& Y_;
        PtrListDictionary<volScalarField>& massTransferRates_;
        PtrListDictionary<fvScalarMatrix>& implicitSources_;
        const List<bool>& active_;
        const volScalarField& alphaRho_;
        const surfaceScalarField& alphaRhoPhi_;
        const fvModels& models_;
        const fvConstraints& constraints_;

    public:
        integrator
        (
            const fvMesh& mesh,
            PtrList<volScalarField>& Y,
            PtrListDictionary<volScalarField>& massTransferRates,
            PtrListDictionary<fvScalarMatrix>& implicitSources,
            const List<bool>& active,
            const word& alphaRhoName,
            const word& alphaRhoPhiName
        );

        virtual ~integrator();

        virtual void update();
        virtual void solve();
        virtual void postUpdate();
    };

    autoPtr<integrator> integratorPtr_;



public:

    // Constructor
    multicomponentBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentBlastThermo
    (
        const speciesTable& species,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );


    //- Destructor
    virtual ~multicomponentBlastThermo();


    virtual const word& phaseName() const
    {
        return phaseName_;
    }

    //- ODE functions

        //- Initialize models
        virtual void initializeModels();

        //- Pre update step
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addDelta(const word& name, tmp<volScalarField>& delta);

        //- Add an external change to a thermodynamic field if present
        //  i.e. mass transfer
        virtual void addSource(const word& name, tmp<fvScalarMatrix>& source);


    // Access

        // Specie properties

            //- Molecular weight of the given specie [kg/kmol]
            virtual scalar Wi(const label speciei) const = 0;

            //- Enthalpy of formation [J/kg]
            virtual scalar Hf(const label speciei) const = 0;

            //- Flame temperature [K]
            virtual scalar flameT(const label speciei) const = 0;


        // Specie thermo properties

            //- Density [kg/m^3]
            virtual scalar rho
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual scalar Cp
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> Cp
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual scalar HE
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> HE
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<volScalarField> HE
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual scalar Hs
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<scalarField> Hs
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Sensible enthalpy [J/kg]
            virtual tmp<volScalarField> Hs
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual scalar Ha
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Enthalpy/Internal energy [J/kg]
            virtual tmp<scalarField> Ha
            (
                const label speciei,
                const scalarField& p,
                const scalarField& T
            ) const = 0;

            //- Absolute enthalpy [J/kg]
            virtual tmp<volScalarField> Ha
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;


        // Specie transport properties

            //- Dynamic viscosity [kg/m/s]
            virtual scalar mu
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Dynamic viscosity [kg/m/s]
            virtual tmp<volScalarField> mu
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual scalar kappa
            (
                const label speciei,
                const scalar p,
                const scalar T
            ) const = 0;

            //- Thermal conductivity [W/m/K]
            virtual tmp<volScalarField> kappa
            (
                const label speciei,
                const volScalarField& p,
                const volScalarField& T
            ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentBlastThermoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
