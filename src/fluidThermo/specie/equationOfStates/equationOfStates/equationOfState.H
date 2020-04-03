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
    Foam::equationOfState

Description
    Simple class to transfer the correct quantities for a temperature based
    equation of state

SourceFiles
    equationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef equationOfState_H
#define equationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class equationOfState Declaration
\*---------------------------------------------------------------------------*/

template<class EoS>
class equationOfState
:
    public EoS
{

public:

    // Constructors

        equationOfState(const dictionary& dict);

    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return word(EoS::typeName());
        }

        //- Is temperature used
        bool temperatureBased() const
        {
            return true;
        }

        //- Return pressure
        scalar p(const scalar& rho, const scalar& e, const scalar& T) const;

        //- Return the speed of sound
        scalar speedOfSound
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const scalar& cv
        ) const;

        //- Return the MieGruniesen coefficient
        scalar Gamma
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const scalar& Cv
        ) const;

        //- Return derivative of pressure w.r.t. specific volume
        scalar dpdv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. internal energy
        scalar dpde
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOfStateI.H"

#ifdef NoRepository
    #include "equationOfState.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
