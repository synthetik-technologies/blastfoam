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
    Foam::constTransport

Description
    transport model based with constant properties

SourceFiles
    constTransport.C
\*---------------------------------------------------------------------------*/

#ifndef constTransportBlast_H
#define constTransportBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class constTransport;

template<class Thermo>
inline constTransport<Thermo> operator+
(
    const constTransport<Thermo>&,
    const constTransport<Thermo>&
);

template<class Thermo>
inline constTransport<Thermo> operator*
(
    const scalar,
    const constTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const constTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class const Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constTransport
:
    public Thermo
{
    // Private data

        //- Constant dynamic viscosity [Pa s]
        scalar mu_;

        //- Prandtl Number []
        scalar Pr_;


        //- Construct from components
        inline constTransport
        (
            const Thermo& t,
            const scalar mu,
            const scalar Pr
        );


public:

    // Constructors

        //- Construct as named copy
        inline constTransport(const word&, const constTransport&);

        //- Construct from dictionary
        constTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<constTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<constTransport> New(const dictionary& dict);

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "const<"
              + word(Thermo::typeName()) + '>';
        }


    // Member Functions

        //- Dynamic viscosity [kg/ms]
        scalar mu
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Thermal conductivity [W/mK]
        scalar kappa
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        scalar alphah
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;


    // Member Operators

        inline void operator+=(const constTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        inline friend constTransport operator+ <Thermo>
        (
            const constTransport&,
            const constTransport&
        );

        inline friend constTransport operator* <Thermo>
        (
            const scalar,
            const constTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constTransportI.H"

#ifdef NoRepository
    #include "constTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
