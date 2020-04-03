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



public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "const<"
              + word(Thermo::typeName()) + '>';
        }

    // Constructor
    constTransport(const dictionary& dict);

    // Member Functions

        //- Dynamic viscosity [kg/ms]
        scalar mu
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return mu_;
        }

        //- Thermal conductivity [W/mK]
        scalar kappa
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Thermo::Cp(rho, e, T)*mu(rho, e, T)/Pr_;
        }

        //- Thermal diffusivity of enthalpy [kg/ms]
        scalar alphah
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return mu(rho, e, T)/Pr_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "constTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
