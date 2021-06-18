/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
     \\/     M anipulation  | Synthetik Applied Technologies
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
    Foam::CpCvModel

Description
    Basic class to return Cp and Cv


SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef CpCvModel_H
#define CpCvModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{
/*---------------------------------------------------------------------------*\
                           Class CpCvModel Declaration
\*---------------------------------------------------------------------------*/

class CpCvModel
{
public:

        virtual ~CpCvModel()
        {}

    // Member Functions

        //- Sensible thermal internal energy [J/kg]
        virtual scalar Ets
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Absolute thermal internal energy [J/kg]
        virtual scalar Eta
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Heat capacity at constant volume [J/kg/K]
        virtual scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Sensible thermal enthalpy [J/kg]
        virtual scalar Hts
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Absolute thermal enthalpy [J/kg]
        virtual scalar Hta
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;

        //- Heat capacity at constant pressure [J/kg/K]
        virtual scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
