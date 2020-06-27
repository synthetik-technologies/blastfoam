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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return EquationOfState::Gamma(rho, e, T, T);
}

template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return
        EquationOfState::speedOfSound
        (
            p,
            rho,
            e,
            T,
            T // dummy variable, Must use MG
        );
}

template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::Es
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return TTable_.reverseLookupY(T, rho);
}


template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::Ea
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return TTable_.reverseLookupY(T, rho);
}


template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::Cv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return 1.0/max(TTable_.dFdY(rho, e), small);
}


template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::Cp
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    scalar cv(Cv(rho, e, T));
    scalar cSqr
    (
        sqr
        (
            EquationOfState::speedOfSound
            (
                EquationOfState::p(rho, e, T),
                rho,
                e,
                T,
                T
            )
        )
    );
    return cv/(1.0 - sqr(EquationOfState::Gamma(rho, e, T, e))*cv*T/cSqr);
}


template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::CpMCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Cp(rho, e, T) - Cv(rho, e, T);
}

template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::CpByCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return Cp(rho, e, T)/Cv(rho, e, T);
}


template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::TRhoE
(
    const scalar& T0,
    const scalar& rho,
    const scalar& e
) const
{
    return TTable_.lookup(rho, e);
}


template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::initializeEnergy
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    if (EquationOfState::temperatureBased())
    {
        return TTable_.reverseLookupY(T, rho);
    }

    if (rho < small)
    {
        return 0.0;
    }
    scalar Eest = 1000.0;
    scalar Enew = 1000.0;
    scalar Etol = Eest*tolerance_;
    int    iter = 0;
    do
    {
        Eest = Enew;

        Enew -=
            (EquationOfState::p(rho, Eest, T) - p)
           /stabilise(EquationOfState::dpde(rho, Eest, T), small);

        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }

    } while (mag(Enew - Eest)/Eest > Etol);

    return Enew;
}

template<class EquationOfState>
Foam::scalar Foam::tabulatedThermo<EquationOfState>::initializeRho
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    NotImplemented;
    return rho;
}
// ************************************************************************* //
