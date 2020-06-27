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
    Foam::eConst

Description
    Thermodynamics model based with constant properties

SourceFiles
    eConstI.H
    eConst.C
\*---------------------------------------------------------------------------*/

#ifndef eConstBlast_H
#define eConstBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class eConst Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class eConst
:
    public EquationOfState
{
    // Private data

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv_;

        //- Heat of formation [J/kg]
        scalar Hf_;



public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "eConst<"
              + word(EquationOfState::typeName()) + '>';
        }

    // Constructor
    eConst(const dictionary& dict);


    // Member Functions

        //- Heat capacity at constant volume of mixture [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Cv_ + EquationOfState::Cv(rho, e, T);
        }

        //- Sensible internal energy [J/kg]
        scalar Es
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Cv(rho, e, T)*T + EquationOfState::E(rho, e, T);
        }

        //- Absolute internal energy [J/kg]
        scalar Ea
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Es(rho, e, T) + Hf_;
        }

        #include "EtoHFuncs.H"

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
            return Cp(rho, e, T)*log(T/Tstd) + EquationOfState::S(p, rho, T);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "eConst.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
