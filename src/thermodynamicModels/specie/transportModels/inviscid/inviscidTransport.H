/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::inviscidTransport

Description
    Transport package with 0 properties

Usage

SourceFiles
    inviscidTransportI.H
    inviscidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef inviscidTransport_H
#define inviscidTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class inviscidTransport;

template<class Thermo>
inline inviscidTransport<Thermo> operator+
(
    const inviscidTransport<Thermo>&,
    const inviscidTransport<Thermo>&
);

template<class Thermo>
inline inviscidTransport<Thermo> operator*
(
    const scalar,
    const inviscidTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const inviscidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                       Class inviscidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class inviscidTransport
:
    public Thermo
{
    // Private Member Functions

        //- Construct from components
        inline inviscidTransport(const Thermo& t);


public:

    // Constructors

        //- Construct as named copy
        inline inviscidTransport(const word&, const inviscidTransport&);

        //- Construct from dictionary
        inviscidTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<inviscidTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<inviscidTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "inviscid<" + Thermo::typeName() + '>';
        }

        //- Return if inviscid
        inline bool inviscid() const
        {
            return true;
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const inviscidTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend inviscidTransport operator+ <Thermo>
        (
            const inviscidTransport&,
            const inviscidTransport&
        );

        friend inviscidTransport operator* <Thermo>
        (
            const scalar,
            const inviscidTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const inviscidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "inviscidTransportI.H"

#ifdef NoRepository
    #include "inviscidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
