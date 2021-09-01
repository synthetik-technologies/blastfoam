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
    Foam::multicomponentSolidBlastThermo

Description
    Multicomponent fluid model with species transport

SourceFiles
    multicomponentSolidBlastThermo.C
    multicomponentSolidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentSolidBlastThermo_H
#define multicomponentSolidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentBlastThermo.H"
#include "mixtureBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class multicomponentBlastSolidMixture Declaration
\*---------------------------------------------------------------------------*/

template<class BasicMixture, class ThermoType>
class multicomponentBlastSolidMixture
:
    public BasicMixture
{
    const speciesMixtureField<ThermoType>& mixtureField_;
    mutable UautoPtr<const ThermoType> mixturePtr_;

public:

    multicomponentBlastSolidMixture
    (
        const speciesMixtureField<ThermoType>& mixture
    )
    :
        mixtureField_(mixture),
        mixturePtr_(nullptr)
    {}

    virtual ~multicomponentBlastSolidMixture()
    {}

    //- Set the local reaction progress variable
    void setCell(const label celli) const
    {
        mixturePtr_.reset(&mixtureField_[celli]);
    }
    void setPatchFace(const label patchi, const label facei) const
    {
        mixturePtr_.reset(&mixtureField_.boundary(patchi, facei));
    }

    //- Calculate internal energy
    virtual scalar HE(const scalar rho, const scalar e, const scalar T) const
    {
        return mixturePtr_->Es(rho, e, T);
    }

    //- Calculate temperature
    virtual scalar TRhoE(const scalar T, const scalar rho, const scalar e) const
    {
        return mixturePtr_->TRhoE(T, rho, e);
    }

    //- Calculate specific heat at constant pressure
    virtual scalar Cp(const scalar rho, const scalar e, const scalar T) const
    {
        return mixturePtr_->Cp(rho, e, T);
    }

    //- Calculate specific heat at constant volume
    virtual scalar Cv(const scalar rho, const scalar e, const scalar T) const
    {
        return mixturePtr_->Cv(rho, e, T);
    }

    //- Calculate thermal diffusivity
    virtual scalar kappa(const scalar rho, const scalar e, const scalar T) const
    {
        return mixturePtr_->kappa(rho, e, T);
    }
};

/*---------------------------------------------------------------------------*\
                Class multicomponentSolidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class multicomponentSolidBlastThermo
:
    public Thermo
{
    //- Mixture
    mutable multicomponentBlastSolidMixture
    <
        typename Thermo::BasicMixture,
        typename Thermo::thermoType
    > localMixture_;

    //- Return a cell mixture
    virtual const typename Thermo::BasicMixture& cellThermoMixture
    (
        const label celli
    ) const
    {
        localMixture_.setCell(celli);
        return localMixture_;
    }

    //- Return a cell mixture
    virtual const typename Thermo::BasicMixture& patchFaceThermoMixture
    (
        const label patchi,
        const label facei
    ) const
    {
        localMixture_.setPatchFace(patchi, facei);
        return localMixture_;
    }


public:

    //- Runtime type information
    TypeName("multicomponent");

    // Constructor
    multicomponentSolidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentSolidBlastThermo
    (
        const HashPtrTable<Thermo, word, string::hash>& thermoData,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~multicomponentSolidBlastThermo();


    // Member Functions

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Calculate internal energy
        virtual tmp<volScalarField> calce() const;

        //- Pre update step
        virtual void update()
        {}

        //- Solve sub-step stepi
        virtual void solve()
        {
            multicomponentBlastThermo::solve();
        }

        //- Post update
        virtual void postUpdate()
        {}

    //- Thermodynamic and transport functions

        //- Calculate density
        virtual void updateRho();

        //- Thermal conductivity [W/m/K]
        virtual tmp<volScalarField> kappa() const;

        //- Thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<vectorField> Kappa
        (
            const label patchi
        ) const;

        //- Dynamic viscosity [kg/m/s]
        virtual tmp<volScalarField> mu
        (
            const label speciei,
            const volScalarField& p,
            const volScalarField& T
        ) const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return Thermo::thermoType::isotropic;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentSolidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
