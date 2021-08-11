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
    Foam::multicomponentFluidBlastThermo

Description
    Multicomponent fluid model with species transport

SourceFiles
    multicomponentFluidBlastThermo.C
    multicomponentFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef multicomponentFluidBlastThermo_H
#define multicomponentFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentBlastThermo.H"
#include "mixtureBlastThermo.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class multicomponentFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class multicomponentFluidBlastThermo
:
    public Thermo
{
public:

    //- Runtime type information
    TypeName("multicomponent");

    // Constructor
    multicomponentFluidBlastThermo
    (
        const word& phaseName,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const word& masterName = word::null
    );

    // Constructor
    multicomponentFluidBlastThermo
    (
        const HashPtrTable<Thermo, word, string::hash>& thermoData,
        const word& phaseName,
        volScalarField& rho,
        volScalarField& e,
        volScalarField& T,
        const dictionary& dict,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~multicomponentFluidBlastThermo();


    // Member Functions

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return the speed of sound
        virtual tmp<volScalarField> speedOfSound() const;

        //- Return the speed of sound for patchi
        virtual tmp<scalarField> speedOfSound(const label patchi) const;

        //- Gamma = dpde []
        virtual tmp<volScalarField> Gamma() const;

        //- Gamma = dpde for patch []
        virtual tmp<scalarField> Gamma(const label patchi) const;

        //- Gamma = dpde for cell []
        virtual scalar Gammai(const label celli) const;


    //- Thermodynamic and transport functions

        //- Calculate density
        virtual void updateRho(const volScalarField& p);

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> pRhoT() const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce(const volScalarField& p) const;

    // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual tmp<volScalarField> mu() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addMulticomponentFluidThermo(Transport, Thermo, EoS)                \
    typedef multicomponentFluidBlastThermo<Transport##Thermo##EoS##blastSpecie>  \
        multicomponent##Transport##Thermo##EoS;                             \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        multicomponent##Transport##Thermo##EoS,                             \
        (Transport##Thermo##EoS##blastSpecie::typeName()).c_str(),          \
        0                                                                   \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        fluidBlastThermo,                                                   \
        multicomponent##Transport##Thermo##EoS,                             \
        multicomponent                                                      \
    );

#define addMulticomponentFluidThermos(Transport, EoS)                       \
    addMulticomponentFluidThermo(Transport, eConst, EoS);                   \
    addMulticomponentFluidThermo(Transport, hConst, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "multicomponentFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
