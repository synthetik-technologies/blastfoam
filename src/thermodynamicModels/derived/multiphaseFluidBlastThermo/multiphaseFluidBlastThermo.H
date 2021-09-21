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
    Foam::multiphaseFluidTherm

Description
    Class to calculate mixture properties of a collection of more than two
    equation of states.

    References:
    \verbatim
        Zheng, H.W., Shu, C., Chew, Y.T., Qin, N.  (2011).
        A solution adaptive simulation of compressible multi-fluid flows
        with general equation of state.
        International Journal for Numerical Methods in Fluids, 67, 616-637.
    \endverbatim

SourceFiles
    multiphaseFluidBlastThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multiphaseFluidBlastThermo_H
#define multiphaseFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseFluidBlastThermo.H"
#include "fluidBlastThermo.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class multiphaseFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

class multiphaseFluidBlastThermo
:
    public fluidBlastThermo
{
// Private data

        //- Names of phases
        wordList phases_;

        //- Volume fraction of phases
        PtrList<volScalarField> volumeFractions_;

        //- Volume weighted density of phases
        PtrList<volScalarField> rhos_;

        //- Equation of state for phases
        PtrList<phaseFluidBlastThermo> thermos_;

        //- Volume weighted density of phases
        PtrList<volScalarField> alphaRhos_;

        //- Residual volume fraction
        dimensionedScalar residualAlpha_;

        //- Pointer to total volume fraction field
        UautoPtr<const volScalarField> sumVfPtr_;


    // Private member functions

        //- Correct thermodynamic fields
        virtual void calculate();

        //- Return the field normalized by the sumVfPtr field if it is set
        tmp<volScalarField> normalise(tmp<volScalarField>& tf) const
        {
            return
                sumVfPtr_.valid()
              ? tf/max(sumVfPtr_(), residualAlpha_)
              : tf;
        }

        //- Return the field normalized by the sumVfPtr field if it is set
        tmp<scalarField> normalise(tmp<scalarField>& tf, const label patchi) const
        {
            return
                sumVfPtr_.valid()
              ? tf/max(sumVfPtr_->boundaryField()[patchi], residualAlpha_.value())
              : tf;
        }

        //- Return the field normalized by the sumVfPtr field if it is set
        tmp<scalarField> normalise
        (
            tmp<scalarField>& tf,
            const labelList& cells
        ) const
        {
            return
                sumVfPtr_.valid()
              ? tf
                /max
                (
                    UIndirectList<scalar>(sumVfPtr_(), cells)(),
                    residualAlpha_.value()
                )
              : tf;
        }

        //- Return a cell normalized by the sumVfPtr field if it is set
        scalar normalise(const scalar f, const label celli) const
        {
            return
                sumVfPtr_.valid()
              ? f/max(sumVfPtr_()[celli], residualAlpha_.value())
              : f;
        }

        //- Return the field normalized by the sumVfPtr field if it is set
        void normalise(volScalarField& f) const
        {
            if (sumVfPtr_.valid())
            {
                f /= max(sumVfPtr_(), residualAlpha_);
            }
        }


public:

    TypeName("multiphaseFluid");


    // Constructor
    multiphaseFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~multiphaseFluidBlastThermo();


     // Member Functions

        //- Correct thermodynamic fields
        virtual void correct();

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
            return "multiphase";
        }

    // Access functions

        //- Non-const access to volume fractions
        PtrList<volScalarField>& volumeFractions()
        {
            return volumeFractions_;
        }

        //- Const access to volume fractions
        const PtrList<volScalarField>& volumeFractions() const
        {
            return volumeFractions_;
        }

        //- Set the pointer to the total volume fraction field
        void setTotalVolumeFractionPtr(const volScalarField& vf);

        //- Non-const access to densities
        PtrList<volScalarField>& rhos()
        {
            return rhos_;
        }

        //- Const access to densities
        const PtrList<volScalarField>& rhos() const
        {
            return rhos_;
        }

        //- Non-const access to thermo model phasei
        virtual phaseFluidBlastThermo& thermo(const label phasei)
        {
            return thermos_[phasei];
        }

        //- Const access to thermo model phasei
        virtual const phaseFluidBlastThermo& thermo(const label phasei) const
        {
            return thermos_[phasei];
        }


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Recompute the density
            virtual void updateRho(const volScalarField&);

            //- Calculate pressure for cell [Pa]
            virtual scalar cellpRhoT(const label celli) const;

            //- Gamma = dpde for cell []
            virtual scalar cellGamma(const label celli) const;

            //- Derivative of pressure w.r.t density for cel
            virtual scalar celldpdRho(const label celli) const;

            //- Derivative of pressure w.r.t internal energy for cell
            virtual scalar celldpde(const label celli) const;


        //- Thermodynamic and transport functions

            using fluidBlastThermo::he;

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

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> Cp
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for a cell [J/kg/K]
            virtual scalar cellCp(const scalar, const label) const;

            //- Heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> Cv
            (
                const scalarField& T,
                const label patchi
            ) const;

            //- Heat capacity at constant volume for a cell [J/kg/K]
            virtual scalar cellCv(const scalar, const label) const;

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
