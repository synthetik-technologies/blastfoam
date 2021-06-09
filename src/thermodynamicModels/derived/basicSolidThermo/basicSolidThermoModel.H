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
    Foam::basicSolidThermoModel

Description
    Non-reacting, single phase solid thermodynamic model

SourceFiles
    basicSolidThermoModel.C
    basicSolidThermoModels.C

\*---------------------------------------------------------------------------*/

#ifndef basicSolidThermoModel_H
#define basicSolidThermoModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solidThermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class basicSolidThermoModel Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class basicSolidThermoModel
:
    public Thermo
{
public:

    //- Runtime type information
    TypeName("basicSolidThermo");

    // Constructor
    basicSolidThermoModel
    (
        const word& phaseName,
        const fvMesh& mesh,
        const dictionary& dict,
        const bool master,
        const word& masterName = word::null
    );

    //- Destructor
    virtual ~basicSolidThermoModel();


    // Member Functions

        //- Pre update
        virtual void update()
        {}

        //- Solve sub-step stepi
        virtual void solve()
        {}

        //- Remove stored fields
        virtual void clearODEFields()
        {}

        //- Post update
        virtual void postUpdate()
        {}

        //- Correct fields
        virtual void correct();


    //- Thermodynamic and transport functions

        //- Calculate density
        virtual void updateRho()
        {}

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
            return Thermo::thermoType::isotropic();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define addSolidThermo(Transport, Thermo, EoS)                              \
    typedef basicSolidThermoModel                                           \
        <                                                                   \
            eThermoModel                                                    \
            <                                                               \
                solidThermoModel,                                           \
                Transport##Thermo##EoS##blastSpecie                         \
            >                                                               \
        > basic##Transport##Thermo##EoS;                                    \
    defineTemplateTypeNameAndDebugWithName                                  \
    (                                                                       \
        basic##Transport##Thermo##EoS,                                      \
        (Transport##Thermo##EoS##blastSpecie::typeName()).c_str(),          \
        0                                                                   \
    );                                                                      \
    addToRunTimeSelectionTable                                              \
    (                                                                       \
        solidThermoModel,                                                   \
        basic##Transport##Thermo##EoS,                                      \
        basicSolid                                                          \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "basicSolidThermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
