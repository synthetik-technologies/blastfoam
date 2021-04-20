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
    Foam::MGEquationOfState

Description

    Simple class to transfer the correct quantities for MieGruniesen based
    equation of state

SourceFiles
    MGEquationOfState.C

\*---------------------------------------------------------------------------*/

#ifndef MGEquationOfState_H
#define MGEquationOfState_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EoS> class MGEquationOfState;

template<class EoS>
inline MGEquationOfState<EoS> operator+
(
    const MGEquationOfState<EoS>&,
    const MGEquationOfState<EoS>&
);

template<class EoS>
inline MGEquationOfState<EoS> operator*
(
    const scalar,
    const MGEquationOfState<EoS>&
);

template<class EoS>
inline MGEquationOfState<EoS> operator==
(
    const MGEquationOfState<EoS>&,
    const MGEquationOfState<EoS>&
);

template<class EoS>
Ostream& operator<<
(
    Ostream&,
    const MGEquationOfState<EoS>&
);


/*---------------------------------------------------------------------------*\
                           Class MGEquationOfState Declaration
\*---------------------------------------------------------------------------*/

template<class EoS>
class MGEquationOfState
:
    public EoS
{

public:

    //- Construct from components
        inline MGEquationOfState
        (
            const EoS& sp
        );

        //- Construct from dictionary
        MGEquationOfState(const dictionary& dict);

        //- Construct as named copy
        inline MGEquationOfState
        (
            const word& name,
            const MGEquationOfState&
        );

        //- Construct and return a clone
        inline autoPtr<MGEquationOfState> clone() const;

        // Selector from dictionary
        inline static autoPtr<MGEquationOfState> New
        (
            const dictionary& dict
        );


    //- Return the instantiated type name
        static word typeName()
        {
            return word(EoS::typeName());
        }


    // Member Functions

        //- Is temperature used
        bool temperatureBased() const
        {
            return false;
        }

        //- Return pressure
        scalar p
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const bool limit = true
        ) const;

        //- Return the speed of sound
        scalar cSqr
        (
            const scalar& p,
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const scalar& Cv
        ) const;

        //- Return derivative of pressure w.r.t. specific volume
        scalar dpdv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. Temperature
        scalar dpdT
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return derivative of pressure w.r.t. internal energy
        scalar dpde
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Return the MieGruniesen coefficient
        scalar Gamma
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T,
            const scalar& Cv
        ) const;

        //- Internal energy correction [J/kg]
        scalar E
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Enthalpy correction [J/kg]
        scalar H
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant volume [J/kg/K]
        scalar Cv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity at constant pressure [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

        //- Heat capacity difference [J/kg/K]
        scalar CpMCv
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const;

    // Member operators

        inline void operator+=(const MGEquationOfState&);
        inline void operator*=(const scalar);


    // Friend operators

        friend MGEquationOfState operator+ <EoS>
        (
            const MGEquationOfState&,
            const MGEquationOfState&
        );

        friend MGEquationOfState operator* <EoS>
        (
            const scalar s,
            const MGEquationOfState&
        );

        friend MGEquationOfState operator== <EoS>
        (
            const MGEquationOfState&,
            const MGEquationOfState&
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MGEquationOfStateI.H"

#ifdef NoRepository
    #include "MGEquationOfState.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
