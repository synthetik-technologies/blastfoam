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
    Foam::keyesTransport

Description
    Transport package using Sutherland's formula for viscosity:

    \verbatim
        mu = mu0*sqrt(T)/(1 + (B/T)*10^(-C/T))
    \endverbatim

    and the thermal conductivity (for gases) is obtained from:

    \verbatim
        kappa = mu*Cp/Pr
    \endverbatim

Usage
    \table
        Property     | Description
        mu0          | Reference viscosity
        B            | Reference temperature [K]
        C            | Power [K]
    \endtable

    Example specification of keyesTransport for air:
    \verbatim
    transport
    {
        mu0             1.458e-06;
        B               110.4;
    }
    \endverbatim

SourceFiles
    keyesTransportI.H
    keyesTransport.C

\*---------------------------------------------------------------------------*/

#ifndef keyesTransport_H
#define keyesTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class keyesTransport;

template<class Thermo>
inline keyesTransport<Thermo> operator+
(
    const keyesTransport<Thermo>&,
    const keyesTransport<Thermo>&
);

template<class Thermo>
inline keyesTransport<Thermo> operator*
(
    const scalar,
    const keyesTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const keyesTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                     Class keyesTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class keyesTransport
:
    public Thermo
{
    // Private Data

        // Sutherland's coefficients
        scalar mu0_, B_, C_;

        // Optional Prandtl number
        scalar rPr_;


    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);

        //- Read coefficient from dictionary
        void readCoeffs(const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        inline keyesTransport
        (
            const Thermo& t,
            const scalar mu0,
            const scalar T0,
            const scalar C,
            const scalar rPr
        );

        //- Construct as named copy
        inline keyesTransport(const word&, const keyesTransport&);

        //- Construct from dictionary
        keyesTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        keyesTransport(const Thermo& t,const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<keyesTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<keyesTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "keyes<" + Thermo::typeName() + '>';
        }

        //- Return if inviscid
        inline bool inviscid() const
        {
            return false;
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

        inline void operator+=(const keyesTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend keyesTransport operator+ <Thermo>
        (
            const keyesTransport&,
            const keyesTransport&
        );

        friend keyesTransport operator* <Thermo>
        (
            const scalar,
            const keyesTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const keyesTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyesTransportI.H"

#ifdef NoRepository
    #include "keyesTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
