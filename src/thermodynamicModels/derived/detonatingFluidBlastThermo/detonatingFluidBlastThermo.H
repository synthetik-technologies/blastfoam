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
    Foam::detonatingFluidBlastThermo

Description
    Detonating fluid model with options to include activation and afterburn
    models

SourceFiles
    detonatingFluidBlastThermo.C
    detonatingFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef detonatingFluidBlastThermo_H
#define detonatingFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"
#include "afterburnModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class detonatingBlastFluidMixture Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMixture, class ThermoType1, class ThermoType2>
class detonatingBlastFluidMixture
:
    public BasicMixture
{
    const ThermoType1& thermo1_;
    const ThermoType2& thermo2_;
    const activationModel& activation_;

    mutable scalar xi_;

public:

    detonatingBlastFluidMixture
    (
        const ThermoType1& thermo1,
        const ThermoType2& thermo2,
        const activationModel& activation
    )
    :
        thermo1_(thermo1),
        thermo2_(thermo2),
        activation_(activation)
    {}

    virtual ~detonatingBlastFluidMixture()
    {}

    //- Set the local reaction progress variable
    void setCell(const label celli) const
    {
        xi_ = activation_.cellLambdaPow(celli);
    }
    void setPatchFace(const label patchi, const label facei) const
    {
        xi_ = activation_.lambdaPow(patchi, facei);
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

    //- Calculate dpde
    virtual scalar Gamma(const scalar rho, const scalar e, const scalar T) const;

    //- Calculate speed of sound squared
    virtual scalar cSqr(const scalar, const scalar, const scalar, const scalar) const;
};

/*---------------------------------------------------------------------------*\
                           Class detonatingFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class detonatingFluidBlastThermo
:
    public Thermo
{
// Private data

    //- Activation model
    autoPtr<activationModel> activation_;

    //- Afterburn model
    autoPtr<afterburnModel> afterburn_;

    //- Mixture
    mutable detonatingBlastFluidMixture
    <
        typename Thermo::BasicMixture,
        typename Thermo::thermoType1,
        typename Thermo::thermoType2
    > mixture_;

    //- Return a cell mixture
    virtual const typename Thermo::BasicMixture& cellThermoMixture
    (
        const label celli
    ) const
    {
        mixture_.setCell(celli);
        return mixture_;
    }

    //- Return a cell mixture
    virtual const typename Thermo::BasicMixture& patchFaceThermoMixture
    (
        const label patchi,
        const label facei
    ) const
    {
        mixture_.setPatchFace(patchi, facei);
        return mixture_;
    }

    //- Return the blending field for celli
    virtual scalar cellx(const label celli) const
    {
        return
            activation_.valid()
          ? activation_->cellLambdaPow(celli)
          : 0.0;
    }

    //- Return the blending field for patchi
    virtual tmp<scalarField> x(const label patchi) const
    {
        return
            activation_.valid()
          ? activation_->lambdaPow(patchi)
          : tmp<scalarField>(new scalarField(this->rho_.size(), 0.0));
    }


public:

    TypeName("detonating");

    // Constructor
    detonatingFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~detonatingFluidBlastThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Post update
        virtual void postUpdate();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;


    //- Thermodynamic and transport functions

        //- Calculate density
        virtual void updateRho(const volScalarField& p);

        //- Calculate pressure for cell [Pa]
        virtual scalar cellpRhoT(const label celli) const;

        //- Gamma = dpde for cell []
        virtual scalar cellGamma(const label celli) const;

        //- Derivative of pressure w.r.t density for cell
        virtual scalar celldpdRho(const label celli) const;

        //- Derivative of pressure w.r.t internal energy for cell
        virtual scalar celldpde(const label celli) const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce(const volScalarField& p) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "detonatingFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
