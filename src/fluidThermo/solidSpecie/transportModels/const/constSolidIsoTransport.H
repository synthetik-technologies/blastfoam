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
    Foam::constSolidIsoTransport

Description
    transport model based with constant properties

SourceFiles
    constSolidIsoTransport.C
\*---------------------------------------------------------------------------*/

#ifndef constSolidIsoTransportBlast_H
#define constSolidIsoTransportBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class constSolidIsoTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class constSolidIsoTransport
:
    public Thermo
{
    // Private data

        //- Constant isotropic thermal conductivity
        scalar kappa_;


public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "constIso<"
              + word(Thermo::typeName()) + '>';
        }

    // Constructor
    constSolidIsoTransport(const dictionary& dict);

    // Member Functions

        //- Is the thermal conductivity isotropic
        bool isotropic() const
        {
            return true;
        }

        //- Dynamic viscosity [kg/ms]
        scalar mu
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return 0.0;
        }

        //- Thermal conductivity [W/mK]
        scalar kappa
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return kappa_;
        }

        //- Thermal conductivity [W/mK]
        vector Kappa
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return vector(kappa_, kappa_, kappa_);
        }

        //- Thermal diffusivity of enthalpy [kg/ms]
        scalar alphah
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return kappa(rho, e, T)/this->Cp(rho, e, T);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "constSolidIsoTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
