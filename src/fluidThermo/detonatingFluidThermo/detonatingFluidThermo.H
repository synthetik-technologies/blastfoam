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
    Foam::detonatingFluidThermo

Description
    Templated class to allow for a detonating material. Material will
    transition from the unreacted state to reacted state using a activation
    model and add the detonation energy.

SourceFiles
    detonatingFluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef detonatingFluidThermo_H
#define detonatingFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidThermoModel.H"
#include "activationModel.H"
#include "afterburnModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class detonatingFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class uThermo, class rThermo>
class detonatingFluidThermo
:
    public fluidThermoModel,
    public uThermo,
    public rThermo
{
// Private data

    //- Activation model
    autoPtr<activationModel> activation_;

    //- Afterburn model
    autoPtr<afterburnModel> afterburn_;


protected:

    //- Protected functions

        //- Return a volScalarField of the given property
        template<class uMethod, class rMethod, class ... Args>
        tmp<volScalarField> volScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            uMethod upsiMethod,
            rMethod rpsiMethod,
            const Args& ... args
        ) const;

        //- Return a volScalarField of the given property for a single state
        template<class Method, class ... Args>
        tmp<volScalarField> singleVolScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Method psiMethod,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a patch
        template<class uMethod, class rMethod, class ... Args>
        tmp<scalarField> patchFieldProperty
        (
            uMethod upsiMethod,
            rMethod rpsiMethod,
            const label patchi,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a cell set
        template<class uMethod, class rMethod, class ... Args>
        tmp<scalarField> cellSetProperty
        (
            uMethod upsiMethod,
            rMethod rpsiMethod,
            const labelList& cells,
            const Args& ... args
        ) const;


public:

    TypeName("detonatingFluidThermo");

    // Constructor
    detonatingFluidThermo
    (
        const word& phaseName,
        volScalarField& p,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const bool master
    );


    //- Destructor
    virtual ~detonatingFluidThermo();


    // Member Functions

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        );

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const labelList& oldIs,
            const label& nOld,
            const labelList& deltaIs,
            const label nDelta
        );

        //- Remove stored fields
        virtual void clearODEFields();

        //- Correct fields
        virtual void correct();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const;


    //- Thermodynamic and transport functions

//         //- Calculate density
//         virtual tmp<volScalarField> calcRho
//         (
//             const volScalarField& p,
//             const volScalarField& e,
//             const volScalarField& T
//         ) const;

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> calcP() const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce() const;

        //- Calculate temperature
        virtual tmp<volScalarField> calcT() const;

        //- Calculate internal energy
        virtual tmp<volScalarField> E() const;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> E
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const label patchi
        ) const;

        //- Calculate internal energy for patchi
        virtual tmp<scalarField> E
        (
            const scalarField& rho,
            const scalarField& e,
            const scalarField& T,
            const labelList& faceCells
        ) const;

        //- Return molecular weight
        virtual tmp<volScalarField> W() const;

        //- Return molecular weight
        virtual scalar Wi(const label celli) const;

        //- Mie Gruniesen coefficient
        virtual tmp<volScalarField> Gamma() const;

        //- Mie Gruniesen coefficient for a patch
        virtual tmp<scalarField> Gamma(const label patchi) const;

        //- Heat capacity at constant pressure of mixture [J/kg/K]
        virtual tmp<volScalarField> Cp() const;

        //- Heat capacity at constant pressure of mixture for a patch [J/kg/K]
        virtual tmp<scalarField> Cp(const label patchi) const;

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

        //- Specific heat ratio for a patch
        virtual tmp<scalarField> CpByCv(const label patchi) const;

        //- Read residual values
        virtual void read(const dictionary& dict)
        {
            residualRho_.read(dict);
            residualAlpha_.read(dict);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addDetonatingFluidThermo(uTransport, uThermo, uEoSType, uEoS, rTransport, rThermo, rEoSType, rEoS) \
    typedef detonatingFluidThermo                                               \
        <                                                                       \
            uTransport##uThermo##uEoSType##uEoS,                                \
            rTransport##rThermo##rEoSType##rEoS                                 \
        >                                                                       \
        detonating##uTransport##uThermo##uEoSType##uEoS##rTransport##rThermo##rEoSType##rEoS; \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        detonating##uTransport##uThermo##uEoSType##uEoS##rTransport##rThermo##rEoSType##rEoS, \
        (                                                                       \
            "detonating<"                                                       \
          + word(uTransport##uThermo##uEoSType##uEoS::typeName())               \
          + ','                                                                 \
          + word(rTransport##rThermo##rEoSType##rEoS::typeName())               \
          + '>'                                                                 \
        ).c_str(),                                                              \
        0                                                                       \
    );                                                                          \
    addToRunTimeSelectionTable                                                  \
    (                                                                           \
        fluidThermoModel,                                                       \
        detonating##uTransport##uThermo##uEoSType##uEoS##rTransport##rThermo##rEoSType##rEoS,                                       \
        detonating                                                               \
    );

#define addDetonatingFluidThermos(uTransport, uEoSType, uEoS, rTransport, rEoSType, rEoS) \
    addDetonatingFluidThermo                                                    \
    (                                                                           \
        uTransport, eConst, uEoSType, uEoS,                                     \
        rTransport, eConst, rEoSType, rEoS                                      \
    );                                                                          \
    addDetonatingFluidThermo                                                    \
    (                                                                           \
        uTransport, hConst, uEoSType, uEoS,                                     \
        rTransport, hConst, rEoSType, rEoS                                      \
    );                                                                          \
    addDetonatingFluidThermo                                                    \
    (                                                                           \
        uTransport, hConst, uEoSType, uEoS,                                     \
        rTransport, eConst, rEoSType, rEoS                                      \
    );                                                                          \
    addDetonatingFluidThermo                                                    \
    (                                                                           \
        uTransport, eConst, uEoSType, uEoS,                                     \
        rTransport, hConst, rEoSType, rEoS                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "detonatingFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
