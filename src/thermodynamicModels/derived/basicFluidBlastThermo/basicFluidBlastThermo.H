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
    Foam::basicFluidBlastThermo

Description
    Non reacting, single phase fluid thermo model

SourceFiles
    basicFluidBlastThermo.C
    basicFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef basicFluidBlastThermo_H
#define basicFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class basicBlastFluidMixture Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMixture, class ThermoType>
class basicBlastFluidMixture
:
    public BasicMixture
{
    const ThermoType& thermo_;

public:

    basicBlastFluidMixture(const ThermoType& thermo)
    :
        thermo_(thermo)
    {}

    virtual ~basicBlastFluidMixture()
    {}

    //- Calculate internal energy
    virtual scalar HE(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.Es(rho, e, T);
    }

    //- Calculate temperature
    virtual scalar TRhoE(const scalar T, const scalar rho, const scalar e) const
    {
        return thermo_.TRhoE(T, rho, e);
    }

    //- Calculate specific heat at constant pressure
    virtual scalar Cp(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.Cp(rho, e, T);
    }

    //- Calculate specific heat at constant volume
    virtual scalar Cv(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.Cv(rho, e, T);
    }

    //- Calculate thermal diffusivity
    virtual scalar kappa(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.kappa(rho, e, T);
    }

    //- Calculate pressure
    virtual scalar pRhoT(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.p(rho, e, T);
    }

    //- Calculate viscosity
    virtual scalar mu(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.mu(rho, e, T);
    }

    //- Calculate dpde
    virtual scalar Gamma(const scalar rho, const scalar e, const scalar T) const
    {
        return thermo_.Gamma(rho, e, T);
    }

    //- Calculate speed of sound squared
    virtual scalar speedOfSound
    (
        const scalar p,
        const scalar rho,
        const scalar e,
        const scalar T
    ) const
    {
        return sqrt(max(thermo_.cSqr(p, rho, e, T), small));
    }

    //- Calculate speed of sound squared
    virtual scalar cSqr
    (
        const scalar p,
        const scalar rho,
        const scalar e,
        const scalar T
    ) const
    {
        return thermo_.cSqr(p, rho, e, T);
    }
};


/*---------------------------------------------------------------------------*\
                           Class basicFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class basicFluidBlastThermo
:
    public Thermo
{
    //- Mixture
    mutable basicBlastFluidMixture
    <
        typename Thermo::BasicMixture,
        typename Thermo::thermoType
    > mixture_;

    //- Return a cell mixture
    virtual const typename Thermo::BasicMixture& cellThermoMixture
    (
        const label celli
    ) const
    {
        return mixture_;
    }

    //- Return a cell mixture
    virtual const typename Thermo::BasicMixture& patchFaceThermoMixture
    (
        const label patchi,
        const label facei
    ) const
    {
        return mixture_;
    }


public:

    //- Runtime type information
    TypeName("basic");

    // Constructor
    basicFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~basicFluidBlastThermo();


    // Member Functions

        //- Pre update
        virtual void update()
        {}

        //- Solve sub-step stepi
        virtual void solve()
        {}

        //- Remove stored fields
        virtual void clearODEFields()
        {}

        //- Post update
        virtual void postUpdate()
        {}

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
    #include "basicFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
