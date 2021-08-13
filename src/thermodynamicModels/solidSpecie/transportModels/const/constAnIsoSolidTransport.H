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
    Foam::constAnIsoSolidTransport

Description
    transport model based with constant properties

SourceFiles
    constAnIsoSolidTransport.C

\*---------------------------------------------------------------------------*/

#ifndef constAnIsoSolidTransportBlast_H
#define constAnIsoSolidTransportBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

template<class Thermo> class constAnIsoSolidTransport;

template<class Thermo>
inline constAnIsoSolidTransport<Thermo> operator+
(
    const constAnIsoSolidTransport<Thermo>&,
    const constAnIsoSolidTransport<Thermo>&
);

template<class Thermo>
inline constAnIsoSolidTransport<Thermo> operator*
(
    const scalar,
    const constAnIsoSolidTransport<Thermo>&
);

/*---------------------------------------------------------------------------*\
                           Class constAnIsoSolidTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constAnIsoSolidTransport
:
    public Thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
        vector kappa_;

        // Private Member Functions

        //- Construct from components
        inline constAnIsoSolidTransport(const Thermo& t, const vector kappa);


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "constAnIso<"
              + word(Thermo::typeName()) + '>';
        }

    // Constructors

        //- Construct as named copy
        inline constAnIsoSolidTransport
        (
            const word&,
            const constAnIsoSolidTransport&
        );

        //- Construct from dictionary
        constAnIsoSolidTransport(const dictionary& dict);


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

        inline void operator+=(const constAnIsoSolidTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend constAnIsoSolidTransport operator+ <Thermo>
        (
            const constAnIsoSolidTransport&,
            const constAnIsoSolidTransport&
        );

        friend constAnIsoSolidTransport operator* <Thermo>
        (
            const scalar,
            const constAnIsoSolidTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "constAnIsoSolidTransportI.H"

#ifdef NoRepository
    #include "constAnIsoSolidTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
