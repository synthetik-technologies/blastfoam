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
    Foam::twoPhaseFluidThermo

Description
    Class to calculate mixture properties of two equation of states

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    twoPhaseFluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseFluidThermo_H
#define twoPhaseFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class twoPhaseFluidThermo Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseFluidThermo
:
    public fluidThermoModel
{
// Private data

        //- Names of phases
        wordList phases_;

        //- Volume fraction
        volScalarField volumeFraction_;

        //- Density of phase 1
        volScalarField rho1_;

        //- Equation of state for phase 1
        autoPtr<fluidThermoModel> thermo1_;

        //- Density of phase 2
        volScalarField rho2_;

        //- Equation of state for phase 2
        autoPtr<fluidThermoModel> thermo2_;


public:


    // Constructor
    twoPhaseFluidThermo
    (
        const word& name,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        volScalarField& p,
        const dictionary& dict,
        const bool master
    );


    //- Destructor
    virtual ~twoPhaseFluidThermo();


    // Member Functions

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const labelList& oldIs,
            const label& nOld,
            const labelList& deltaIs,
            const label nDelta
        );

        //- Remove stored fields
        virtual void clearODEFields();

        //- Correct fields
        virtual void correct();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const;

    // Access functions

        //- Const access to volume fraction
        const volScalarField& volumeFraction() const
        {
            return volumeFraction_;
        }

        //- Non-cont access to volume fraction
        volScalarField& volumeFraction()
        {
            return volumeFraction_;
        }

        //- Const reference to first equation of state
        const fluidThermoModel& thermo1() const
        {
            return thermo1_();
        }

        //- Non-const reference to first equation of state
        fluidThermoModel& thermo1()
        {
            return thermo1_();
        }

        //- Const reference to second equation of state
        const fluidThermoModel& thermo2() const
        {
            return thermo2_();
        }

        //- Non-const reference to second equation of state
        fluidThermoModel& thermo2()
        {
            return thermo2_();
        }


    //- Thermodynamic and transport functions

//         //- Calculate density
//         virtual tmp<volScalarField> calcRho
//         (
//             const volScalarField& p,
//             const volScalarField& e,
//             const volScalarField& T
//         ) const;

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> calcP() const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce() const;

        //- Calculate temperature
        virtual tmp<volScalarField> calcT() const;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> TRhoE
        (
            const scalarField& T,
            const scalarField& e,
            const label patchi
        ) const;

        //- Calculate internal energy for celli
        virtual scalar TRhoEi
        (
            const scalar& T,
            const scalar& e,
            const label celli
        ) const;

        //- Calculate internal energy
        virtual tmp<volScalarField> E() const;

        //- Calculate internal energy for patchi
        virtual tmp<volScalarField> e
        (
            const volScalarField& rho,
            const volScalarField& e,
            const volScalarField& T
        ) const;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> e
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> e
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const labelList& faceCells
        ) const;

        //- Return molecular weight
        virtual tmp<volScalarField> W() const;

        //- Return molecular weight
        virtual tmp<scalarField> W(const label patchi) const;

        //- Return molecular weight
        virtual scalar Wi(const label celli) const;

        //- Mie Gruniesen coefficient
        virtual tmp<volScalarField> Gamma() const;

        //- Mie Gruniesen coefficient for a patch
        virtual tmp<scalarField> Gamma(const label patchi) const;

        //- Heat capacity at constant pressure of mixture [J/kg/K]
        virtual tmp<volScalarField> Cp() const;

        //- Heat capacity at constant pressure of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cp(const label patchi) const;

        //- Heat capacity at constant pressure of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cp
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity at constant pressure of mixture for a cell [J/kg/K]
        virtual scalar Cpi(const label celli) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        virtual tmp<volScalarField> Cv() const;

        //- Heat capacity at constant volume of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cv(const label patchi) const;

        //- Heat capacity at constant volume of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cv
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const;

        //- Heat capacity at constant volume of mixture for a cell [J/kg/K]
        virtual scalar Cvi(const label celli) const;

        //- Specific heat ratio
        virtual tmp<volScalarField> CpByCv() const;

        //- Specific heat ratio for a patch
        virtual tmp<scalarField> CpByCv(const label patchi) const;

        //- Specific heat ratio for patchi
        virtual tmp<scalarField> CpByCv
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
