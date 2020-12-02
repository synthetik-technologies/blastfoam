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
    Foam::detonatingSolidThermo

Description
    Templated class to allow for a detonating material. Material will
    transition from the unreacted state to reacted state using a activation
    model and add the detonation energy.

SourceFiles
    detonatingSolidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef detonatingSolidThermo_H
#define detonatingSolidThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "activationModel.H"
#include "afterburnModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class detonatingSolidThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class detonatingSolidThermo
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

    TypeName("detonatingSolidThermo");

    // Constructor
    detonatingSolidThermo
    (
        const word& phaseName,
        const fvMesh& mesh,
        const dictionary& dict,
        const bool master,
        const word& masterName = word::null
    );

    //- Initialize models
    virtual void initializeModels();


    //- Destructor
    virtual ~detonatingSolidThermo();


    // Member Functions

        //- Pre update
        virtual void update();

        //- Solve sub-step stepi
        virtual void solve();

        //- Remove stored fields
        virtual void clearODEFields();

        //- Post update
        virtual void postUpdate();

        //- Correct fields
        virtual void correct();

        //- Return energy source
        virtual tmp<volScalarField> ESource() const;

        //- Thermal conductivity [W/m/K]
        virtual tmp<volVectorField> Kappa() const;

        //- Anisotropic thermal conductivity [W/m/K]
        virtual tmp<vectorField> Kappa
        (
            const label patchi
        ) const;

        //- Return true if thermal conductivity is isotropic
        virtual bool isotropic() const
        {
            return
                Thermo::thermoType1::isotropic()
             && Thermo::thermoType2::isotropic();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addDetonatingSolidThermo(Transport, Thermo, uEoSType, uEoS, rEoSType, rEoS) \
    typedef detonatingSolidThermo                                           \
        <                                                                   \
            blendedThermoModel                                              \
            <                                                               \
                solidThermoModel,                                           \
                Transport##Thermo##uEoS##blastSpecie,                       \
                Transport##Thermo##rEoS##rBlastSpecie                       \
            >                                                               \
        > detonating##Transport##Thermo##uEoS##rEoS;                        \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        detonating##Transport##Thermo##uEoS##rEoS,                          \
        (                                                                   \
            "detonating<"                                                   \
          + word(Transport##Thermo##uEoS##blastSpecie::typeName())          \
          + ','                                                             \
          + word(Transport##Thermo##rEoS##rBlastSpecie::typeName())         \
          + '>'                                                             \
        ).c_str(),                                                          \
        0                                                                   \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        solidThermoModel,                                                   \
        detonating##Transport##Thermo##uEoS##rEoS,                          \
        detonatingSolid                                                     \
    );

#define addDetonatingSolidThermos(Transport)                                \
    addDetonatingSolidThermo                                                \
    (                                                                       \
        Transport, eConst,                                                  \
        equationOfState, rhoConst,                                          \
        equationOfState, rhoConst                                           \
    );                                                                      \
    addDetonatingSolidThermo                                                \
    (                                                                       \
        Transport, hConst,                                                  \
        equationOfState, rhoConst,                                          \
        equationOfState, rhoConst                                           \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "detonatingSolidThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
