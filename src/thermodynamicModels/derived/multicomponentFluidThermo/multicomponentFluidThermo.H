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
    Multicomponent fluid model with species transport

SourceFiles
    multicomponentFluidThermo.C
    multicomponentFluidThermos.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentFluidThermo_H
#define multicomponentFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentThermoModels.H"
#include "mixtureThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multicomponentFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class multicomponentFluidThermo
:
    public mixtureThermoModel<multicomponentFluidThermoModel, ThermoType>
{
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

    // Constructor
    multicomponentFluidThermo
    (
        const HashPtrTable<ThermoType, word, string::hash>& thermoData,
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

        //- Pre update
        virtual void update()
        {}

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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addMulticomponentFluidThermo(Transport, Thermo, EoSType, EoS)       \
    typedef multicomponentFluidThermo<Transport##Thermo##EoS##blastSpecie>  \
        multicomponent##Transport##Thermo##EoS;                             \
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

#ifdef NoRepository
    #include "multicomponentFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
