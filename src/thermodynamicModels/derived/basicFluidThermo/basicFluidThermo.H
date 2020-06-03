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
    Foam::basicFluidThermo

Description
    Template class to define a runtime selectable fluidThermo model

SourceFiles
    basicFluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef basicFluidThermo_H
#define basicFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class basicFluidThermo
:
    public Thermo
{
public:

    //- Runtime type information
    TypeName("basicFluidThermo");

    // Constructor
    basicFluidThermo
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
    virtual ~basicFluidThermo();


    // Member Functions

        //- Solve sub-step stepi
        virtual void solve
        (
            const label stepi,
            const scalarList& ai,
            const scalarList& bi
        )
        {}

        //- Set old lists and fluxes (initialization of fields)
        virtual void setODEFields
        (
            const label nSteps,
            const labelList& oldIs,
            const label& nOld,
            const labelList& deltaIs,
            const label nDelta
        )
        {}

        //- Remove stored fields
        virtual void clearODEFields()
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

        //- Calculate thermodynamic pressure
        virtual tmp<volScalarField> calcP() const;

        //- Calculate internal energy from pressure and density
        virtual tmp<volScalarField> calce() const;

        //- Read residual values
        virtual void read(const dictionary& dict)
        {
            this->residualRho_.read(dict);
            this->residualAlpha_.read(dict);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addFluidThermo(Transport, Thermo, EoSType, EoS)                         \
    typedef basicFluidThermo<eThermoModel<fluidThermoModel, Transport##Thermo##EoSType##EoS>>     \
        basicFluidThermo##Transport##Thermo##EoSType##EoS;                      \
    defineTemplateTypeNameAndDebugWithName                                      \
    (                                                                           \
        basicFluidThermo##Transport##Thermo##EoSType##EoS,                      \
        (Transport##Thermo##EoSType##EoS::typeName()).c_str(),                  \
        0                                                                       \
    );                                                                          \
    addToRunTimeSelectionTable                                                  \
    (                                                                           \
        fluidThermoModel,                                                       \
        basicFluidThermo##Transport##Thermo##EoSType##EoS,                      \
        basic                                                                   \
    );

#define addFluidThermos(Transport, EoSType, EoS)                            \
    addFluidThermo(Transport, eConst, EoSType, EoS);                        \
    addFluidThermo(Transport, hConst, EoSType, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
