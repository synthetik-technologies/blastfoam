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
    Detonating fluid model with options to include activation and afterburn
    models

SourceFiles
    detonatingFluidThermo.C
    detonatingFluidThermos.C

\*---------------------------------------------------------------------------*/

#ifndef detonatingFluidThermo_H
#define detonatingFluidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"
#include "afterburnModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class detonatingFluidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class detonatingFluidThermo
:
    public Thermo
{
// Private data

    //- Activation model
    autoPtr<activationModel> activation_;

    //- Afterburn model
    autoPtr<afterburnModel> afterburn_;

    //- Return the blending field for celli
    virtual scalar xi(const label celli) const
    {
        return activation_->lambdaPowi(celli);
    }

    //- Return the blending field for patchi
    virtual tmp<scalarField> x(const label patchi) const
    {
        return activation_->lambdaPow(patchi);
    }


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
        const bool master,
        const word& masterName = word::null
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~detonatingFluidThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

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

        //- Calculate density
        virtual void updateRho();

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

#define addDetonatingFluidThermo(Transport, Thermo, uEoS, rEoS)             \
    typedef detonatingFluidThermo                                           \
        <                                                                   \
            blendedThermoModel                                              \
            <                                                               \
                fluidThermoModel,                                           \
                Transport##Thermo##uEoS##blastSpecie,                       \
                Transport##Thermo##rEoS##blastSpecie                        \
            >                                                               \
        > detonating##Transport##Thermo##uEoS##rEoS;                        \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        detonating##Transport##Thermo##uEoS##rEoS,                          \
        (                                                                   \
            "detonating<"                                                   \
          + word(Transport##Thermo##uEoS##blastSpecie::typeName())          \
          + ','                                                             \
          + word(Transport##Thermo##rEoS##blastSpecie::typeName())          \
          + '>'                                                             \
        ).c_str(),                                                          \
        0                                                                   \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        fluidThermoModel,                                                   \
        detonating##Transport##Thermo##uEoS##rEoS,                          \
        detonating                                                          \
    );

#define addDetonatingFluidThermos(Transport, uEoS, rEoS)                    \
    addDetonatingFluidThermo(Transport, eConst, uEoS, rEoS);                \
//     addDetonatingFluidThermo(Transport, hConst, uEoS, rEoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "detonatingFluidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
