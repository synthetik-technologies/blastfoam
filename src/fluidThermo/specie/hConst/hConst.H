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
    Foam::hConst

Description
    Thermodynamics model based with constant properties

SourceFiles
    hConstI.H
    hConst.C
\*---------------------------------------------------------------------------*/

#ifndef hConstBlast_H
#define hConstBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class hConst Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class hConst
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp_;

        //- Heat of formation [J/kg]
        scalar Hf_;



public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "hConst<"
              + word(EquationOfState::typeName()) + '>';
        }

    // Constructor
    hConst(const dictionary& dict);


    // Member Functions

        //- Heat capacity at constant pressure of mixture [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Cp_ + EquationOfState::Cp(rho, e, T);
        }

        //- Sensible enthalpy [J/kg]
        scalar Hs
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Cp(rho, e, T)*T + EquationOfState::H(rho, e, T);
        }

        //- Absolute enthalpy [J/kg]
        scalar Ha
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Hs(rho, e, T) + Hf_;
        }

        #include "HtoEFuncs.H"

        //- Heat of formation [J/kg]
        scalar Hf
        (
            const scalar& rho,
            const scalar& T
        ) const
        {
            return Hf_;
        }

        //- Entropy [J/kg/K]
        scalar S
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Cp(rho, e, T)*log(T/Tstd) + EquationOfState::S(rho, T);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "hConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
