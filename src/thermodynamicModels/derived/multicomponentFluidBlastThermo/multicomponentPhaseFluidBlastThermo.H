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
    Foam::multicomponentPhaseFluidBlastThermo

Description
    Multicomponent phase fluid model with species transport

SourceFiles
    multicomponentPhaseFluidBlastThermo.C
    multicomponentPhaseFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentPhaseFluidBlastThermo_H
#define multicomponentPhaseFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtureBlastThermo.H"
#include "phaseFluidBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class multicomponentPhaseFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class multicomponentPhaseFluidBlastThermo
:
    public Thermo
{

public:

    //- Runtime type information
    TypeName("multicomponent");

    // Constructor
    multicomponentPhaseFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentPhaseFluidBlastThermo
    (
        const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~multicomponentPhaseFluidBlastThermo();


    // Member Functions

        //- Add contribution to mixture temperature
        virtual void calculate
        (
            const volScalarField& alpha,
            const volScalarField& he,
            const volScalarField& T,
            volScalarField& Cp,
            volScalarField& Cv,
            volScalarField& mu,
            volScalarField& alphah,
            volScalarField& pXiSum,
            volScalarField& XiSum
        );

        //- Add contribution to speed of sound
        virtual void calculateSpeedOfSound
        (
            const volScalarField& alpha,
            volScalarField& cSqrRhoXiSum
        );

        //- Pre update step
        virtual void update()
        {
            multicomponentBlastThermo::update();
        }

        //- Solve sub-step stepi
        virtual void solve()
        {
            multicomponentBlastThermo::solve();
        }

        //- Explicit post update
        virtual void solveExplicit()
        {
            multicomponentBlastThermo::solveExplicit();
        }

        //- Store explicit deltas
        virtual void storeExplicit()
        {
            multicomponentBlastThermo::storeExplicit();
        }

        //- Implicit post update
        virtual void solveImplicit()
        {
            multicomponentBlastThermo::solveImplicit();
        }

        //- Clear stored data
        virtual void clear()
        {
            multicomponentBlastThermo::clear();
        }

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return energy source given the current volume fraction
        virtual tmp<volScalarField> ESource
        (
            const volScalarField& alpha
        ) const;

        //- Return initial energy source
        virtual tmp<volScalarField> initESource() const;

        //- Return initial energy source given the current volume fraction
        virtual tmp<volScalarField> initESource
        (
            const volScalarField& alpha
        ) const;


    //- Thermodynamic and transport functions

        //- Return if inviscid
        virtual bool inviscid() const;

        //- Calculate density
        virtual void updateRho
        (
            const volScalarField& alpha,
            const volScalarField& p
        );

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce(const volScalarField& p) const;

        //- Calculate internal energy for a cell
        virtual scalar calcCelle(const scalar p, const label) const;

        //- Calculate pressure for cell [Pa]
        virtual scalar cellpRhoT
        (
            const label celli,
            const bool limit = true
        ) const;

        //- Calculate pressure for cell [Pa]
        virtual scalar patchFacepRhoT
        (
            const label patchi,
            const label facei,
            const bool limit = true
        ) const;

        //- Calculate Gamma = dpde
        virtual tmp<volScalarField> Gamma() const;

        //- Gamma = dpde for cell []
        virtual scalar cellGamma(const label celli) const;

        //- Gamma = dpde for cell []
        virtual scalar patchFaceGamma(const label patchi, const label facei) const;

        //- Calculate the speed of sound
        virtual tmp<volScalarField> calcSpeedOfSound() const;

        //- Calculate the speed of sound in a cell
        virtual scalar cellSpeedOfSound
        (
            const scalar p,
            const label celli
        ) const;

        //- Derivative of pressure w.r.t density for cell
        virtual scalar celldpdRho(const label celli) const;

        //- Derivative of pressure w.r.t internal energy for cell
        virtual scalar celldpde(const label celli) const;

        //- Derivative of pressure w.r.t temperature for cell
        virtual scalar celldpdT(const label celli) const;

        //- Is the equation of state temperature based
        virtual bool temperatureBased() const
        {
            return Thermo::thermoType::temperatureBased();
        }

        //- basicThermo stubs
        virtual IOdictionary& properties() { NotImplemented; return const_cast<IOdictionary&>(reinterpret_cast<const IOdictionary&>(dictionary::null)); }
        virtual void correct() { NotImplemented; }
        virtual bool incompressible() const { NotImplemented; return false; }
        virtual bool isochoric() const { NotImplemented; return false; }
        virtual Switch dpdt() const { NotImplemented; return Switch(false); }

        virtual const volScalarField& T() const { return this->T_; }
        virtual volScalarField& T() { return const_cast<volScalarField&>(this->T_); }
        virtual const volScalarField& he() const { return this->e_; }
        virtual volScalarField& he() { return const_cast<volScalarField&>(this->e_); }

        virtual tmp<volScalarField> rho() const { return this->rho_; }
        virtual volScalarField& rho() { return this->rho_; }
        virtual tmp<scalarField> rho(const label patchi) const { return this->rho_.boundaryField()[patchi]; }

        //- basicThermo implementation
        virtual const volScalarField& Cp() const { return this->CpRef(); }
        virtual const volScalarField& Cv() const { return this->CvRef(); }
        virtual const volScalarField& Cpv() const { return this->CvRef(); }
        virtual const volScalarField& kappa() const { return this->kappaRef(); }

        // basicSpecieBlastMixture stubs
        virtual scalar cellp(const label celli) const { return this->p()[celli]; }
        virtual scalar pi(const label speciei, const scalar rho, const scalar e, const scalar T) const { return this->speciesData_[speciei].pRhoT(rho, e, T); }
        virtual tmp<volScalarField> pi(const label speciei, const volScalarField& rho, const volScalarField& e, const volScalarField& T) const { return this->volScalarFieldSpecieProperty(speciei, "p", dimPressure, &Thermo::thermoType::pRhoT, rho, e, T); }
        virtual scalar dpdRhoi(const label speciei, const scalar rho, const scalar e, const scalar T) const { return -this->speciesData_[speciei].dpdv(rho, e, T)/sqr(max(rho, SMALL)); }
        virtual scalar dpdTi(const label speciei, const scalar rho, const scalar e, const scalar T) const { return this->speciesData_[speciei].dpdT(rho, e, T); }
        virtual scalar mui(const label speciei, const scalar rho, const scalar e, const scalar T) const { return this->speciesData_[speciei].mu(rho, e, T); }
        virtual tmp<volScalarField> mui(const label speciei, const volScalarField& rho, const volScalarField& e, const volScalarField& T) const { return this->volScalarFieldSpecieProperty(speciei, "mu", dimensionSet(1, -1, -1, 0, 0), &Thermo::thermoType::mu, rho, e, T); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentPhaseFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
