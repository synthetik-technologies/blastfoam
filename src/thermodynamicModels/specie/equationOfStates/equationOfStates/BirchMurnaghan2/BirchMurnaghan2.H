/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
    Foam::BirchMurnaghan2

Description
    Second order BirchMurnaghan equation of state for a solid

    References:
    \verbatim
        Zheng, Z., Zhao, J. (2016).
        Unreacted equation of states of typical energetic material under static
        compression: A review
        Chin. Phys. B, 25 (7), 1-11.
    \endverbatim

SourceFiles
    BirchMurnaghan2.C

\*---------------------------------------------------------------------------*/

#ifndef BirchMurnaghan2Blast_H
#define BirchMurnaghan2Blast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class BirchMurnaghan Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class BirchMurnaghan2
:
    public Specie
{
// Private data

    //- Reference density
    scalar rho0_;

    //- Reference pressure
    scalar pRef_;

    //- Bulk modulus
    scalar K0_;

    //- Derivative of bulk modulus w.r.t. pressure
    scalar K0Prime_;

    //- Mie Grunseisen coefficient
    scalar Gamma_;


public:

    // Constructors

        //- Multiphase constuctor
        BirchMurnaghan2(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "BirchMurnaghan2<"
              + word(Specie::typeName()) + '>';
        }

        //- Is this a solid equation of state
        static bool solid()
        {
            return true;
        }

        //- Return pressure
        scalar p
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Gamma function
        scalar Gamma
        (
            const scalar& rho,
            const scalar& T,
            const scalar& cv
        ) const;

        //- Speed of sound
        scalar speedOfSound
        (
            const scalar& p,
            const scalar& rho,
            const scalar& T,
            const scalar& cv
        ) const;

        //- Return derivative of pressure w.r.t density
        scalar dpdv
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar S
        (
            const scalar& p,
            const scalar& rho,
            const scalar& T
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BirchMurnaghan2I.H"

#ifdef NoRepository
    #include "BirchMurnaghan2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
