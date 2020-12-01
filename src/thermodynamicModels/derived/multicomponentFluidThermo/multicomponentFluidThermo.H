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
    Foam::multicomponentFluidThermo

Description
    Template class to define a runtime selectable fluidThermo model

SourceFiles
    multicomponentFluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentFluidThermo_H
#define multicomponentFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidThermoModel.H"
#include "hashedWordList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multicomponentFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class multicomponentFluidThermo
:
    public fluidThermoModel
{
//- Private data

    //- Name of species
    hashedWordList species_;

    //- Mass fractions
    PtrList<volScalarField> Ys_;

    //- Species data
    PtrList<Thermo> speciesData_;

    //- Stored thermo state at cell centres
    PtrList<Thermo> volMixture_;

    //- Stored thermo state at boundaries
    PtrList<PtrList<Thermo>> faceMixture_;

    //- Temporary storage for the cell/face mixture thermo data
    mutable Thermo mixture_;

    //- Constant pointer to flux
    const surfaceScalarField* alphaRhoPhiPtr_;

    //- Constant point to alphaRho field
    const volScalarField* alphaRhoPtr_;

    //- ODE fields
    PtrList<PtrList<volScalarField>> YsOld_;
    PtrList<PtrList<volScalarField>> deltaAlphaRhoYs_;

    // Private Member Functions

        //- Construct the species data from the given dictionary and return the
        //  data for the first specie to initialise the mixture thermo data
        const Thermo& constructSpeciesData(const dictionary& thermoDict);

        const Thermo& cellMixture(const label celli) const;

        const Thermo& patchFaceMixture
        (
            const label patchi,
            const label facei
        ) const;

        void updateMixture();


public:

    //- Runtime type information
    TypeName("multicomponentFluidThermo");

    // Constructor
    multicomponentFluidThermo
    (
        const word& phaseName,
        volScalarField& p,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const bool master,
        const word& masterName = word::null
    );


    //- Destructor
    virtual ~multicomponentFluidThermo();


    // Member Functions

        //- Initialize models
        virtual void initializeModels();

        //- Pre update
        virtual void update()
        {}

        //- Solve sub-step stepi
        virtual void solve();

        //- Remove stored fields
        virtual void clearODEFields();

        //- Post update
        virtual void postUpdate()
        {}

        //- Correct fields
        virtual void correct();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const;


    //- Thermodynamic and transport functions

        //- Calculate thermodynamic pressure for a patch
        virtual tmp<scalarField> calcP(const label patchi) const;

        //- Calculate thermodynamic pressure
        virtual scalar calcPi(const label celli) const;

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

        //- Return molecular weight for patchi
        virtual tmp<scalarField> W(const label patchi) const;

        //- Return molecular weight
        virtual scalar Wi(const label celli) const;

        //- Mie Gruniesen coefficient
        virtual tmp<volScalarField> Gamma() const;

        //- Mie Gruniesen coefficient for a patch
        virtual tmp<scalarField> Gamma(const label patchi) const;

        //- Mie Gruniesen coefficient for a cell
        virtual scalar Gammai(const label celli) const;

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

        //- Specific heat ratio for patchi
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

#define addMulticomponentFluidThermo(Transport, Thermo, EoSType, EoS)           \
    typedef multicomponentFluidThermo< Transport##Thermo##EoSType##EoS##blastSpecie> \
        multicomponentFluid##Transport##Thermo##EoSType##EoS;                   \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        multicomponentFluid##Transport##Thermo##EoSType##EoS,                   \
        (Transport##Thermo##EoSType##EoS##blastSpecie::typeName()).c_str(),     \
        0                                                                       \
    );                                                                          \
    addToRunTimeSelectionTable                                                  \
    (                                                                           \
        fluidThermoModel,                                                       \
        multicomponentFluid##Transport##Thermo##EoSType##EoS,                   \
        multicomponent                                                          \
    );

#define addMulticomponentFluidThermos(Transport, EoSType, EoS)                  \
    addMulticomponentFluidThermo(Transport, eConst, EoSType, EoS);              \
    addMulticomponentFluidThermo(Transport, hConst, EoSType, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
