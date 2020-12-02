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
    Foam::mixtureThermoModel

Description
    Template class to define a runtime selectable fluidThermo model

SourceFiles
    mixtureThermoModel.C

\*---------------------------------------------------------------------------*/

#ifndef mixtureThermoModel_H
#define mixtureThermoModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtureThermoModels.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class mixtureThermoModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class ThermoType>
class mixtureThermoModel
:
    public BasicThermo
{
protected:
//- Protected data

    //- Species data
    PtrList<Thermo> speciesData_;

    //- Temporary storage for the cell/face mixture thermo data
    mutable Thermo mixture_;

    //- Stored thermo state at cell centres
    PtrList<Thermo> cellMixtures_;

    //- Stored thermo state at boundaries
    PtrList<PtrList<Thermo>> faceMixtures_;


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

        //- Return a scalarField of the given property on a patch
        template<class Method, class ... Args>
        tmp<scalarField> patchFieldProperty
        (
            Method psiMethod,
            const label patchi,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a cell set
        template<class Method, class ... Args>
        tmp<scalarField> cellSetProperty
        (
            Method psiMethod,
            const labelList& cells,
            const Args& ... args
        ) const;


public:

    //- Runtime type information
    TypeName("mixtureThermoModel");

    // Constructor
    mixtureThermoModel
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
    virtual ~mixtureThermoModel();


    // Member Functions

        //- Thermodynamic and transport functions

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define definemixtureThermoModel(Transport, Thermo, EoSType, EoS)    \
    typedef mixtureThermoModel< Transport##Thermo##EoS##blastSpecie> \
        multicomponent##Transport##Thermo##EoS;

#define definemixtureThermoModels(Transport, EoSType, EoS)           \
    definemixtureThermoModel(Transport, eConst, EoSType, EoS);       \
    definemixtureThermoModel(Transport, hConst, EoSType, EoS);


#define addmixtureThermoModel(Transport, Thermo, EoSType, EoS)       \
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

#define addmixtureThermoModels(Transport, EoSType, EoS)              \
    addmixtureThermoModel(Transport, eConst, EoSType, EoS);          \
    addmixtureThermoModel(Transport, hConst, EoSType, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "mixtureThermoModelI.H"

#ifdef NoRepository
    #include "mixtureThermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
