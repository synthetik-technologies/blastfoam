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
    Foam::constIsoSolidTransport

Description
    transport model based with constant properties

SourceFiles
    constIsoSolidTransport.C
\*---------------------------------------------------------------------------*/

#ifndef constIsoSolidTransportBlast_H
#define constIsoSolidTransportBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

template<class Thermo> class constIsoSolidTransport;

template<class Thermo>
inline constIsoSolidTransport<Thermo> operator+
(
    const constIsoSolidTransport<Thermo>&,
    const constIsoSolidTransport<Thermo>&
);


template<class Thermo>
inline constIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const constIsoSolidTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class constIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constIsoSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
        scalar kappa_;

    // Private Member Functions

        //- Construct from components
        inline constIsoSolidTransport(const Thermo& t, const scalar kappa);

public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "constIso<"
              + word(Thermo::typeName()) + '>';
        }

    // Constructors

        //- Construct as named copy
        inline constIsoSolidTransport
        (
            const word&,
            const constIsoSolidTransport&
        );

        //- Construct from dictionary
        constIsoSolidTransport(const dictionary& dict);

    // Member Functions

        //- Is the thermal conductivity isotropic
        static const bool isotropic = true;

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

        //- Thermal conductivity [W/mK]
        vector Kappa
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

        inline void operator+=(const constIsoSolidTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend constIsoSolidTransport operator+ <Thermo>
        (
            const constIsoSolidTransport&,
            const constIsoSolidTransport&
        );

        friend constIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const constIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constIsoSolidTransportI.H"

#ifdef NoRepository
    #include "constIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
