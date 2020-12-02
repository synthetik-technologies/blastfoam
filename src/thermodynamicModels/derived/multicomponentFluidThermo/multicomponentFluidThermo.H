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

#include "mixtureThermoModels.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multicomponentFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class multicomponentFluidThermo
:
    public mixtureFluidThermoModel
{
//- Private data

    typedef multicomponentFluidThermo<Thermo> reactionThermo;

    //- Species data
    PtrList<Thermo> speciesData_;

    //- Temporary storage for the cell/face mixture thermo data
    mutable Thermo mixture_;

    //- Chemistry reader
    autoPtr<chemistryReader<Thermo>> chemReader_;

    //- Reactions
    PtrList<Reaction<Thermo>> reactions_;

    //- Table of species composition
    speciesCompositionTable speciesComposition_;

    //- Chemisty solver
    autoPtr<ode<StandardChemistryModel<reactionThermo, Thermo>>> chemistry_;

    //- Stored thermo state at cell centres
    PtrList<Thermo> volMixture_;

    //- Stored thermo state at boundaries
    PtrList<PtrList<Thermo>> faceMixture_;


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

        //- Return a volScalarField of the given property
        template<class Method, class ... Args>
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const Args& ... args
        ) const;


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
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Remove stored fields
        virtual void clearODEFields();

        //- Post update
        virtual void postUpdate();

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


    //- Return

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return Thermo::typeName();
        }

        //- Return the table of thermos
        inline const PtrList<Thermo>& speciesData() const;

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }

        //- Return the table of thermos
        inline const PtrList<Reaction<Thermo>>& reactions() const;

        //- Return the table of species
        inline const speciesTable& species() const;

        //- Does the mixture include this specie?
        inline bool contains(const word& specieName) const;

        //- Return true for active species
        inline bool active(label speciei) const;

        //- Return the bool list of active species
        inline const List<bool>& active() const;

        //- Set speciei active
        inline void setActive(label speciei);

        //- Set speciei inactive
        inline void setInactive(label speciei);

        //- Return the mass-fraction fields
        inline PtrList<volScalarField>& Y();

        //- Return the const mass-fraction fields
        inline const PtrList<volScalarField>& Y() const;

        //- Return the mass-fraction field for a specie given by index
        inline volScalarField& Y(const label i);

        //- Return the const mass-fraction field for a specie given by index
        inline const volScalarField& Y(const label i) const;

        //- Return the mass-fraction field for a specie given by name
        inline volScalarField& Y(const word& specieName);

        //- Return the const mass-fraction field for a specie given by name
        inline const volScalarField& Y(const word& specieName) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineMulticomponentFluidThermo(Transport, Thermo, EoSType, EoS)    \
    typedef multicomponentFluidThermo< Transport##Thermo##EoS##blastSpecie> \
        multicomponent##Transport##Thermo##EoS;

#define defineMulticomponentFluidThermos(Transport, EoSType, EoS)           \
    defineMulticomponentFluidThermo(Transport, eConst, EoSType, EoS);       \
    defineMulticomponentFluidThermo(Transport, hConst, EoSType, EoS);


#define addMulticomponentFluidThermo(Transport, Thermo, EoSType, EoS)       \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        multicomponent##Transport##Thermo##EoS,                             \
        (Transport##Thermo##EoS##blastSpecie::typeName()).c_str(),          \
        0                                                                   \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        fluidThermoModel,                                                   \
        multicomponent##Transport##Thermo##EoS,                             \
        multicomponent                                                      \
    );

#define addMulticomponentFluidThermos(Transport, EoSType, EoS)              \
    addMulticomponentFluidThermo(Transport, eConst, EoSType, EoS);          \
    addMulticomponentFluidThermo(Transport, hConst, EoSType, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentFluidThermoI.H"

#ifdef NoRepository
    #include "multicomponentFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
