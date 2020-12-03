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
    Foam::reactingFluidThermo

Description
    Template class to define a runtime selectable fluidThermo model

SourceFiles
    reactingFluidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef reactingFluidThermo_H
#define reactingFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentFluidThermo.H"
#include "chemistryReader.H"
#include "BasicChemistryModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reactingFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class reactingFluidThermo
:
    public speciesTable,
    public autoPtr<chemistryReader<ThermoType>>,
    public multicomponentFluidThermo<ThermoType>,
    public PtrList<Reaction<ThermoType>>
{
// Private data

    //- Table of species composition
    speciesCompositionTable speciesComposition_;


    //- Pointer to chemistry model
    autoPtr<BasicChemistryModel<multicomponentFluidThermoModel>> chemistryPtr_;


public:

    //- Runtime type information
    TypeName("reactingFluidThermo");

    // Constructor
    reactingFluidThermo
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
    virtual ~reactingFluidThermo();


    // Member Functions

        //- Post update
        virtual void postUpdate();

        Reaction<ThermoType>& operator[](const label i)
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        const Reaction<ThermoType>& operator[](const label i) const
        {
            return PtrList<Reaction<ThermoType>>::operator[](i);
        }

        //- Table of species composition
        const speciesCompositionTable& specieComposition() const
        {
            return speciesComposition_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addreactingFluidThermo(Transport, Thermo, EoSType, EoS)       \
    typedef reactingFluidThermo<Transport##Thermo##EoS##blastSpecie>  \
        reacting##Transport##Thermo##EoS;                             \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        reacting##Transport##Thermo##EoS,                             \
        (Transport##Thermo##EoS##blastSpecie::typeName()).c_str(),          \
        0                                                                   \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        fluidThermoModel,                                                   \
        reacting##Transport##Thermo##EoS,                             \
        reacting                                                      \
    );

#define addreactingFluidThermos(Transport, EoSType, EoS)              \
    addreactingFluidThermo(Transport, eConst, EoSType, EoS);          \
    addreactingFluidThermo(Transport, hConst, EoSType, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reactingFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
