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
    Foam::reactingPhaseFluidBlastThermo

Description
    Reacting multicomponent phase fluid model with species transport

SourceFiles
    reactingPhaseFluidBlastThermo.C
    reactingPhaseFluidBlastThermos.C

\*---------------------------------------------------------------------------*/

#ifndef reactingPhaseFluidBlastThermo_H
#define reactingPhaseFluidBlastThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "multicomponentPhaseFluidBlastThermo.H"
#include "blastCombustionModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class reactingPhaseFluidBlastThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class reactingPhaseFluidBlastThermo
:
    public multicomponentPhaseFluidBlastThermo<Thermo>
{
// Private member data

    //- Combustion model
    autoPtr<blastCombustionModel> combustion_;

    //- Solve combustion at all ode steps
    bool odeCombustion_;

    //- Heat release rate [J/m^3/s]
    volScalarField Qdot_;


public:

    //- Runtime type information
    TypeName("reacting");

    // Constructor
    reactingPhaseFluidBlastThermo
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    // Constructor
    reactingPhaseFluidBlastThermo
    (
        const HashPtrTable<typename Thermo::thermoType, word, string::hash>& thermoData,
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~reactingPhaseFluidBlastThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Return energy source given the current volume fraction
        virtual tmp<volScalarField> ESource
        (
            const volScalarField& alpha
        ) const;

        //- Return temperature [K]
        virtual const volScalarField& T() const { return multicomponentPhaseFluidBlastThermo<Thermo>::T(); }

        //- Return density [kg/m^3]
        virtual tmp<volScalarField> rho() const { return multicomponentPhaseFluidBlastThermo<Thermo>::rho(); }

        //- Return internal energy/enthalpy [J/kg]
        virtual const volScalarField& he() const { return multicomponentPhaseFluidBlastThermo<Thermo>::he(); }

        //- Return heat capacity at constant pressure [J/kg/K]
        virtual const volScalarField& CpRef() const { return multicomponentPhaseFluidBlastThermo<Thermo>::CpRef(); }

        //- Return heat capacity at constant volume [J/kg/K]
        virtual const volScalarField& CvRef() const { return multicomponentPhaseFluidBlastThermo<Thermo>::CvRef(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "reactingPhaseFluidBlastThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
