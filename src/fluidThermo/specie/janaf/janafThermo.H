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
    Foam::janaf

Description
    JANAF tables based thermodynamics package templated
    into the equation of state.

SourceFiles
    janafThermo.C
\*---------------------------------------------------------------------------*/

#ifndef janafBlast_H
#define janafBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalar.H"
#include "FixedList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class janaf Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class janaf
:
    public EquationOfState
{
public:

    // Public data

        static const int nCoeffs_ = 7;
        typedef FixedList<scalar, nCoeffs_> coeffArray;


private:

    // Private Data

        // Temperature limits of applicability of functions
        scalar Tlow_, Thigh_, Tcommon_;

        coeffArray highCpCoeffs_;
        coeffArray lowCpCoeffs_;


    // Private Member Functions

        //- Return the coefficients corresponding to the given temperature
        const coeffArray& coeffs(const scalar& T) const
        {
            if (T < Tcommon_)
            {
                return lowCpCoeffs_;
            }
            else
            {
                return highCpCoeffs_;
            }
        }



public:

    //- Return the instantiated type name
        static word typeName()
        {
            return
                "janaf<"
              + word(EquationOfState::typeName()) + '>';
        }

    // Constructor
    janaf(const dictionary& dict);


    // Member Functions

        //- Heat capacity at constant pressure of mixture [J/kg/K]
        scalar Cp
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            const coeffArray& a = coeffs(T);
            return
                ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0])
              + EquationOfState::Cp(rho, e, T);
        }

        //- Sensible enthalpy [J/kg]
        scalar Hs
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            return Ha(rho, e, T) - Hf(rho, e, T);
        }

        //- Absolute enthalpy [J/kg]
        scalar Ha
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            const coeffArray& a = coeffs(T);
            return
            (
                ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
              + a[5]
            ) + EquationOfState::H(rho, e, T);
        }

        #include "HtoEFuncs.H"

        //- Heat of formation [J/kg]
        scalar Hf
        (
            const scalar& rho,
            const scalar& e,
            const scalar& T
        ) const
        {
            const coeffArray& a = lowCpCoeffs_;
            return
            (
                (
                    (
                        ((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd
                      + a[1]/2.0
                    )*Tstd
                  + a[0]
                )*Tstd + a[5]
            );
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
            const coeffArray& a = coeffs(T);
            return
            (
                (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*log(T)
              + a[6]
            ) + EquationOfState::S(p, rho, e, T);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "janafThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
