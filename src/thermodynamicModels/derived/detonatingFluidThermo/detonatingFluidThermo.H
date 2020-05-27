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

#define addDetonatingFluidThermo(uTransport, uThermo, uEoSType, uEoS, rTransport, rThermo, rEoSType, rEoS) \
    typedef detonatingFluidThermo                                               \
        <                                                                       \
            blendedThermoModel                                                  \
            <                                                                   \
                fluidThermoModel,                                               \
                uTransport##uThermo##uEoSType##uEoS,                            \
                rTransport##rThermo##rEoSType##rEoS                             \
            >                                                                   \
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
