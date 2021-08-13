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
    Foam::singlePhaseSolidBlastThermo

Description
    Class to calculate mixture properties of two equation of states

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-Solid flows
        with general equation of state.
        International Journal for Numerical Methods in Solids, 67, 616-637.
    \endverbatim

SourceFiles
    singlePhaseSolidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseSolidBlastThermo_H
#define singlePhaseSolidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidBlastThermo.H"
#include "phaseSolidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class singlePhaseSolidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseSolidBlastThermo
:
    public solidBlastThermo
{
// Private data

        //- Equation of state for phase 1
        autoPtr<phaseSolidBlastThermo> thermoPtr_;


public:

    TypeName("singlePhaseSolid");

    // Constructor
    singlePhaseSolidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );

    // Constructor
    singlePhaseSolidBlastThermo
    (
        const word& thermoType,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~singlePhaseSolidBlastThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Correct fields
        virtual void correct();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return residual values
        virtual dimensionedScalar residualAlpha() const
        {
            return thermoPtr_->residualAlpha();
        }
        virtual dimensionedScalar residualRho() const
        {
            return thermoPtr_->residualRho();
        }

    // Access functions

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return thermoPtr_->isotropic();
        }

        //- Const reference to first equation of state
        const basicBlastThermo& thermo(const label i) const
        {
            return thermoPtr_();
        }

        //- Non-const reference to first equation of state
        basicBlastThermo& thermo(const label i)
        {
            return thermoPtr_();
        }


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Recompute the density
            virtual void updateRho();

            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const;

            //- Anisotropic thermal conductivity for patch [W/m/K]
            virtual tmp<vectorField> Kappa(const label) const;


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

        //- Read residual values
        virtual void read(const dictionary& dict)
        {
            thermoPtr_->read(dict);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
