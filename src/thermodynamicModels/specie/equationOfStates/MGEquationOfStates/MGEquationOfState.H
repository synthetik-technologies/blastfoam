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
    Foam::MGEquationOfState

Description

    Simple class to transfer the correct quantities for MieGruniesen based
    equation of state

SourceFiles
    MGEquationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef MGEquationOfState_H
#define MGEquationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class MGEquationOfState Declaration
\*---------------------------------------------------------------------------*/

class MGEquationOfState
{
protected:

    //- Return the MieGruniesen coefficient
    virtual scalar GammaMG
    (
        const scalar rho,
        const scalar e
    ) const = 0;

    //- Return the deviation from ideal pressure
    virtual scalar Pi
    (
        const scalar rho,
        const scalar e
    ) const = 0;

    //- Return the nonideal speed of sound contribution
    virtual scalar delta
    (
        const scalar p,
        const scalar rho,
        const scalar e
    ) const = 0;


public:

        virtual ~MGEquationOfState()
        {}

    // Member Functions

        //- Is temperature used
        bool temperatureBased() const
        {
            return false;
        }

        //- Return pressure
        virtual scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const bool limit = true
        ) const;

        //- Return the speed of sound
        virtual scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar Cv
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return the MieGruniesen coefficient
        scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar Cv
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        virtual scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        virtual scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Heat capacity difference [J/kg/K]
//         virtual scalar CpMCv
//         (
//             const scalar& rho,
//             const scalar& e,
//             const scalar& T
//         ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MGEquationOfStateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
