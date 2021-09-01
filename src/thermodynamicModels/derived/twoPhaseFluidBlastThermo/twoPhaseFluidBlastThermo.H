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
    Foam::twoPhaseFluidBlastThermo

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
    twoPhaseFluidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef twoPhaseFluidBlastThermo_H
#define twoPhaseFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidBlastThermo.H"
#include "phaseFluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class twoPhaseBlastFluidMixture Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseBlastFluidMixture
:
    public blastFluidMixture
{
    const volScalarField& alpha_;
    const phaseFluidBlastThermo& thermo1_;
    const phaseFluidBlastThermo& thermo2_;
    const volScalarField& rho1_;
    const volScalarField& rho2_;

    mutable scalar alpha1i_;
    mutable scalar alpha2i_;
    mutable UautoPtr<const blastPhaseFluidMixture> mixture1Ptr_;
    mutable UautoPtr<const blastPhaseFluidMixture> mixture2Ptr_;
    mutable scalar rho1i_;
    mutable scalar rho2i_;
public:

    twoPhaseBlastFluidMixture
    (
        const volScalarField& alpha,
        const phaseFluidBlastThermo& thermo1,
        const phaseFluidBlastThermo& thermo2
    )
    :
        alpha_(alpha),
        thermo1_(thermo1),
        thermo2_(thermo2),
        rho1_(thermo1_.rho()),
        rho2_(thermo2_.rho())
    {}

    virtual ~twoPhaseBlastFluidMixture()
    {}

    //- Set the local reaction progress variable
    void setCell(const label celli) const
    {
        alpha1i_ = alpha_[celli];
        alpha2i_ = 1.0 - alpha1i_;
        rho1i_ = rho1_[celli];
        rho2i_ = rho2_[celli];
        mixture1Ptr_.reset(&thermo1_.cellThermoMixture(celli));
        mixture2Ptr_.reset(&thermo2_.cellThermoMixture(celli));
    }
    void setPatchFace(const label patchi, const label facei) const
    {
        alpha1i_ = alpha_.boundaryField()[patchi][facei];
        alpha2i_ = 1.0 - alpha1i_;
        rho1i_ = rho1_.boundaryField()[patchi][facei];
        rho2i_ = rho2_.boundaryField()[patchi][facei];
        mixture1Ptr_.reset(&thermo1_.patchFaceThermoMixture(patchi, facei));
        mixture2Ptr_.reset(&thermo2_.patchFaceThermoMixture(patchi, facei));
    }

    //- Calculate internal energy
    virtual scalar HE(const scalar, const scalar, const scalar) const;

    //- Calculate temperature
    virtual scalar TRhoE(const scalar, const scalar, const scalar) const;

    //- Calculate specific heat at constant pressure
    virtual scalar Cp(const scalar, const scalar, const scalar) const;

    //- Calculate specific heat at constant volume
    virtual scalar Cv(const scalar, const scalar, const scalar) const;

    //- Calculate thermal diffusivity
    virtual scalar kappa(const scalar, const scalar, const scalar) const;

    //- Calculate pressure
    virtual scalar pRhoT(const scalar, const scalar, const scalar) const;

    //- Calculate viscosity
    virtual scalar mu(const scalar, const scalar, const scalar) const;

    //- Calculate speed of sound
    virtual scalar speedOfSound(const scalar, const scalar, const scalar, const scalar) const;

    //- Calculate speed of sound squared
    virtual scalar cSqr(const scalar, const scalar, const scalar, const scalar) const;
};


/*---------------------------------------------------------------------------*\
                           Class twoPhaseFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class twoPhaseFluidBlastThermo
:
    public fluidBlastThermo
{
// Private data

        //- Names of phases
        wordList phases_;

        //- Volume fraction
        volScalarField volumeFraction_;

        //- Density of phase 1
        volScalarField rho1_;

        //- Equation of state for phase 1
        autoPtr<phaseFluidBlastThermo> thermo1_;

        //- Density of phase 2
        volScalarField rho2_;

        //- Equation of state for phase 2
        autoPtr<phaseFluidBlastThermo> thermo2_;

        //- Mixture
        mutable twoPhaseBlastFluidMixture mixture_;

    //- Return a cell mixture
    virtual const blastFluidMixture& cellThermoMixture
    (
        const label celli
    ) const
    {
        mixture_.setCell(celli);
        return mixture_;
    }

    //- Return a cell mixture
    virtual const blastFluidMixture& patchFaceThermoMixture
    (
        const label patchi,
        const label facei
    ) const
    {
        mixture_.setPatchFace(patchi, facei);
        return mixture_;
    }


public:

    TypeName("twoPhaseFluid");

    // Constructor
    twoPhaseFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~twoPhaseFluidBlastThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Calculate internal energy
        virtual tmp<volScalarField> calce(const volScalarField&) const;

        //- Return the thermo name
        virtual word thermoName() const
        {
            return "twoPhase";
        }

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
        const basicBlastThermo& thermo(const label i) const
        {
            return i == 0 ? thermo1_() : thermo2_();
        }

        //- Non-const reference to first equation of state
        basicBlastThermo& thermo(const label i)
        {
            return i == 0 ? thermo1_() : thermo2_();
        }


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Calculate pressure for cell [Pa]
            virtual scalar cellpRhoT(const label celli) const;

            //- Gamma = dpde for cell []
            virtual scalar cellGamma(const label celli) const;

            //- Derivative of pressure w.r.t density for cel
            virtual scalar celldpdRho(const label celli) const;

            //- Derivative of pressure w.r.t internal energy for cell
            virtual scalar celldpde(const label celli) const;

            //- Recompute the density
            virtual void updateRho(const volScalarField&);


        //- Thermodynamic and transport functions

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual const volScalarField& he() const
            {
                return blastThermo::he();
            }

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
            virtual scalar cellTHE
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
            virtual scalar cellCp(const scalar, const label) const;

            //- Heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> Cv() const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for a cell [J/kg/K]
            virtual scalar cellCv(const scalar, const label) const;

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
            virtual scalar cellW(const label celli) const;


        //- Read residual values
        virtual void read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
