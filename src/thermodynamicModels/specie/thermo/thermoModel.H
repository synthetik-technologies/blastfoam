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
    Foam::thermoModel

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

SourceFiles
    thermoI.H
    thermo.C

\*---------------------------------------------------------------------------*/

#ifndef thermoModelBlast_H
#define thermoModelBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class thermoModel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class thermoModel
:
    public ThermoType
{
protected:

    // Protected data

        //- Convergence criteria for temperature
        scalar tolerance_;

        //- Maximum number of iterations
        label maxIter_;


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return  word(ThermoType::typeName());
        }

    // Constructor
    thermoModel(const dictionary& dict);


    // Member Functions

        //- Return Mie Gruniesen coefficient
        scalar Gamma
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return speed of sound
        scalar speedOfSound
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return specific heat ratio
        scalar CpByCv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return temperature
        scalar TRhoE
        (
            const scalar& T0,
            const scalar& rho,
            const scalar& e
        ) const;

        //- Initialize internal energy
        scalar initializeEnergy
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Initialize internal energy
        scalar initializeRho
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineThermoType(Transport, Thermo, EoSType, EoS)                   \
typedef                                                                     \
    Transport                                                               \
    <                                                                       \
        thermoModel                                                         \
        <                                                                   \
            Thermo                                                          \
            <                                                               \
                EoSType                                                     \
                <                                                           \
                    EoS                                                     \
                    <                                                       \
                        blastSpecie                                         \
                    >                                                       \
                >                                                           \
            >                                                               \
        >                                                                   \
    > Transport##Thermo##EoSType##EoS;

#define defineThermoTypes(Transport, EoSType, EoS)                          \
    defineThermoType(Transport, eConst, EoSType, EoS);                      \
    defineThermoType(Transport, hConst, EoSType, EoS);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thermoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
