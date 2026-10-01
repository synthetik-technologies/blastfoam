/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::eJWLC

Description
    Energy-based Mie-Grüneisen JWL equation of state with the full
    three-term (C-form) isentrope pressure and an analytic isentrope
    temperature reference.

    Pressure:
    \f[
        p(\rho, e) =
            A \exp(-R_1 V)
          + B \exp(-R_2 V)
          + C V^{-(\omega+1)}
          + \omega \rho \bigl(e - e_s(V)\bigr)
    \f]

    where \f$ V = \rho_0/\rho \f$ and the isentrope cold-curve energy is:
    \f[
        e_s(V) = v_0 \!\left[
            \frac{A}{R_1} e^{-R_1 V}
          + \frac{B}{R_2} e^{-R_2 V}
          + \frac{C}{\omega} V^{-\omega}
        \right]
    \f]

    The isentrope temperature reference:
    \f[
        T_s(V) = T_\text{CJ} \left(\frac{V_\text{CJ}}{V}\right)^\omega
    \f]

    is encoded in the EOS energy correction so that when paired with
    \c eConst (Tref = 0, Esref = 0, Cv matching this EOS's Cv_) the
    temperature recovers correctly on and off the isentrope:
    \f[
        T(\rho, e) = T_s(V) + \frac{e - e_s(V)}{C_v}
    \f]

    Pair with:
    \verbatim
        thermodynamics
        {
            Cv      <same value as Cv in equationOfState>;
            Hf      0;
            Tref    0;
            Esref   0;
        }
    \endverbatim

SourceFiles
    eJWLC.C

\*---------------------------------------------------------------------------*/

#ifndef eJWLCBlast_H
#define eJWLCBlast_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtr.H"

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class eJWLC;

template<class Specie>
inline eJWLC<Specie> operator+
(
    const eJWLC<Specie>&,
    const eJWLC<Specie>&
);

template<class Specie>
inline eJWLC<Specie> operator*
(
    const scalar,
    const eJWLC<Specie>&
);

template<class Specie>
inline eJWLC<Specie> operator==
(
    const eJWLC<Specie>&,
    const eJWLC<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const eJWLC<Specie>&
);


/*---------------------------------------------------------------------------*\
                           Class eJWLC Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class eJWLC
:
    public Specie
{
// Private data

        //- Reference density [kg/m^3]
        scalar rho0_;

        //- Mie-Grüneisen coefficient [-]
        scalar omega_;

        //- JWL pressure coefficients [Pa]
        scalar A_;
        scalar B_;
        scalar C_;

        //- JWL exponents [-]
        scalar R1_;
        scalar R2_;

        //- Isentrope reference temperature [K] at the CJ state
        scalar T_CJ_;

        //- Reference density at the CJ state [kg/m^3]
        //  Used to compute V_CJ = rho0_ / rho_CJ_
        scalar rho_CJ_;

        //- Reference heat capacity at constant volume [J/kg/K]
        //  Must match the Cv supplied to the paired caloric EOS (eConst)
        scalar Cv_;

        //- Energy offset [J/kg]
        scalar E0_;

    // Private helpers

        //- Dimensionless volume V = rho0/rho
        inline scalar V(const scalar rho) const;

        //- Isentrope pressure p_s(V) [Pa]
        inline scalar ps(const scalar V) const;

        //- Isentrope cold-curve energy e_s(V) [J/kg]
        inline scalar es(const scalar V) const;

        //- Isentrope temperature T_s(V) [K]
        inline scalar Ts(const scalar V) const;


public:

    // Constructors

        //- Construct from components
        inline eJWLC
        (
            const Specie& sp,
            const scalar rho0,
            const scalar omega,
            const scalar A,
            const scalar B,
            const scalar C,
            const scalar R1,
            const scalar R2,
            const scalar T_CJ,
            const scalar rho_CJ,
            const scalar Cv,
            const scalar E0
        );

        //- Construct from dictionary
        eJWLC(const dictionary& dict);

        //- Construct as named copy
        inline eJWLC
        (
            const word& name,
            const eJWLC&
        );

        //- Construct and return a clone
        inline autoPtr<eJWLC> clone() const;

        // Selector from dictionary
        inline static autoPtr<eJWLC> New
        (
            const dictionary& dict
        );


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "eJWLC<"
            + word(Specie::typeName_()) + '>';
        }

        //- Is temperature used as primary variable? No — energy-based
        static bool temperatureBased()
        {
            return false;
        }

        //- Is this a solid equation of state?
        static bool solid()
        {
            return false;
        }

        //- Return pressure [Pa]
        inline scalar p
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const bool limit = true
        ) const;

        //- Return unlimited pressure using (rho, T) — isentrope only
        //  Note: for compatibility; ignores T, uses isentrope
        inline scalar pRhoT
        (
            const scalar rho,
            const scalar T
        ) const;

        //- Return density given temperature and pressure (not implemented)
        inline scalar rho(const scalar p, const scalar T) const;

        //- Return internal energy given pressure and density [J/kg]
        inline scalar ePRho(const scalar p, const scalar rho) const;

        //- Return Mie-Grüneisen coefficient Gamma = omega
        inline scalar Gamma
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return speed of sound squared [m^2/s^2]
        inline scalar cSqr
        (
            const scalar p,
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar cv
        ) const;

        //- Return dp/dv at constant e [Pa*kg/m^3]
        inline scalar dpdv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return dp/de at constant rho [Pa*kg/J]
        inline scalar dpde
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Return dp/dT — zero for energy-based EOS
        inline scalar dpdT
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- EOS energy correction E_eos(rho) [J/kg]
        //  Returns e_s(V) - Cv_ * T_s(V) so that eConst with Tref=0
        //  recovers T = T_s(V) + (e - e_s(V))/Cv on the isentrope
        inline scalar E
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Enthalpy correction [J/kg]
        inline scalar H
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- EOS contribution to Cv — zero (caloric model provides Cv)
        inline scalar Cv
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- EOS contribution to Cp — zero
        inline scalar Cp
        (
            const scalar rho,
            const scalar e,
            const scalar T
        ) const;

        //- Cp - Cv contribution from the EOS [J/kg/K]
        inline scalar CpMCv
        (
            const scalar rho,
            const scalar e,
            const scalar T,
            const scalar CpCv
        ) const;

        //- Entropy contribution integral Cp/T (not implemented)
        inline scalar Sp(const scalar p, const scalar T) const;

        //- Entropy contribution integral Cv/T (not implemented)
        inline scalar Sv(const scalar p, const scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const eJWLC&);
        inline void operator*=(const scalar);


    // Friend operators

        friend eJWLC operator+ <Specie>
        (
            const eJWLC&,
            const eJWLC&
        );

        friend eJWLC operator* <Specie>
        (
            const scalar s,
            const eJWLC&
        );

        friend eJWLC operator== <Specie>
        (
            const eJWLC&,
            const eJWLC&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const eJWLC&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "eJWLCI.H"

#ifdef NoRepository
    #include "eJWLC.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
