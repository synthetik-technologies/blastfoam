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
    Foam::blendedThermoModel

Description
    Templated class to allow for a detonating material. Material will
    transition from the unreacted state to reacted state using a activation
    model and add the detonation energy.

SourceFiles
    blendedThermoModel.C

\*---------------------------------------------------------------------------*/

#ifndef blendedThermoModel_H
#define blendedThermoModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blendedThermoModel Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermo, class Thermo1, class Thermo2>
class blendedThermoModel
:
    public fluidThermoModel,
    public Thermo1,
    public Thermo2
{
protected:

    typedef Thermo1 thermoType1;
    typedef Thermo2 thermoType2;

    // Protected functions

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

        //- Return a volScalarField of the given property
        template<class Method1, class Method2, class ... Args>
        tmp<volScalarField> blendedVolScalarFieldProperty
        (
            const word& psiName,
            const dimensionSet& psiDim,
            Method1 psiMethod1,
            Method2 psiMethod2,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a patch
        template<class Method1, class Method2, class ... Args>
        tmp<scalarField> blendedPatchFieldProperty
        (
            Method1 psiMethod1,
            Method2 psiMethod2,
            const label patchi,
            const Args& ... args
        ) const;

        //- Return a scalarField of the given property on a cell set
        template<class Method1, class Method2, class ... Args>
        tmp<scalarField> blendedCellSetProperty
        (
            Method1 psiMethod1,
            Method2 psiMethod2,
            const labelList& cells,
            const Args& ... args
        ) const;

        //- Return the blending field for celli
        virtual scalar xi(const label celli) const = 0;

        //- Return the blending field for patchi
        virtual tmp<scalarField> x(const label patchi) const = 0;


public:

    // Constructor

        //- Construct with input field
        blendedThermoModel
        (
            const word& phaseName,
            volScalarField& p,
            volScalarField& rho,
            volScalarField& e,
            volScalarField& T,
            const dictionary& dict,
            const dictionary& dict1,
            const dictionary& dict2,
            const bool master = true
        );

        //- Construct and add fields to db
        blendedThermoModel
        (
            const word& phaseName,
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& dict1,
            const dictionary& dict2,
            const bool master = true
        );


    //- Destructor
    virtual ~blendedThermoModel();



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

#define addblendedThermoModel(uTransport, Thermo1, uEoSType, uEoS, rTransport, Thermo2, rEoSType, rEoS) \
    typedef blendedThermoModel                                               \
        <                                                                       \
            uTransport##Thermo1##uEoSType##uEoS,                                \
            rTransport##Thermo2##rEoSType##rEoS                                 \
        >                                                                       \
        detonating##uTransport##Thermo1##uEoSType##uEoS##rTransport##Thermo2##rEoSType##rEoS; \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        detonating##uTransport##Thermo1##uEoSType##uEoS##rTransport##Thermo2##rEoSType##rEoS, \
        (                                                                       \
            "detonating<"                                                       \
          + word(uTransport##Thermo1##uEoSType##uEoS::typeName())               \
          + ','                                                                 \
          + word(rTransport##Thermo2##rEoSType##rEoS::typeName())               \
          + '>'                                                                 \
        ).c_str(),                                                              \
        0                                                                       \
    );                                                                          \
    addToRunTimeSelectionTable                                                  \
    (                                                                           \
        fluidThermoModel,                                                       \
        detonating##uTransport##Thermo1##uEoSType##uEoS##rTransport##Thermo2##rEoSType##rEoS,                                       \
        detonating                                                               \
    );

#define addblendedThermoModels(uTransport, uEoSType, uEoS, rTransport, rEoSType, rEoS) \
    addblendedThermoModel                                                    \
    (                                                                           \
        uTransport, eConst, uEoSType, uEoS,                                     \
        rTransport, eConst, rEoSType, rEoS                                      \
    );                                                                          \
    addblendedThermoModel                                                    \
    (                                                                           \
        uTransport, hConst, uEoSType, uEoS,                                     \
        rTransport, hConst, rEoSType, rEoS                                      \
    );                                                                          \
    addblendedThermoModel                                                    \
    (                                                                           \
        uTransport, hConst, uEoSType, uEoS,                                     \
        rTransport, eConst, rEoSType, rEoS                                      \
    );                                                                          \
    addblendedThermoModel                                                    \
    (                                                                           \
        uTransport, eConst, uEoSType, uEoS,                                     \
        rTransport, hConst, rEoSType, rEoS                                      \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "blendedThermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
