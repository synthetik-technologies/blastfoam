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

#include "thermoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::thermoModel<ThermoType>::thermoModel(const dictionary& dict)
:
    ThermoType(dict),
    tolerance_(dict.lookupOrDefault("tolerance", 1e-6)),
    maxIter_(dict.lookupOrDefault("maxIter", 100))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::Gamma
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return ThermoType::Gamma(rho, e, T, ThermoType::Cv(rho, e, T));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::speedOfSound
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return ThermoType::speedOfSound(p, rho, e, T, ThermoType::Cv(rho, e, T));
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::CpByCv
(
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    return ThermoType::Cp(rho, e, T)/ThermoType::Cv(rho, e, T);
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::TRhoE
(
    const scalar& T0,
    const scalar& rho,
    const scalar& e
) const
{
    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tolerance_;
    int    iter = 0;
    do
    {
        Test = Tnew;
        Tnew =
            Test
          - (ThermoType::Es(rho, e, Test) - e)/ThermoType::Cv(rho, e, Test);
        Tnew = max(Tnew, small);

    } while (mag(Tnew - Test) > Ttol && iter++ < maxIter_);

    return Tnew;
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::initializeEnergy
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    if (ThermoType::temperatureBased())
    {
        scalar E = ThermoType::Es(rho, e, T);
        return E;
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
            (ThermoType::p(rho, Eest, T) - p)
           /stabilise(ThermoType::dpde(rho, Eest, T), small);

        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }

    } while (mag(Enew - Eest)/Eest > Etol);

    return Enew;
}

template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::initializeRho
(
    const scalar& p,
    const scalar& rho,
    const scalar& e,
    const scalar& T
) const
{
    //- Simple method to calculate initial density
    //  Should be modified to solve the 2D problem for density and internal energy
    scalar Rhoest = max(1e-4, rho);
    scalar Rhonew = Rhoest;
    scalar Rhotol = Rhoest*tolerance_;
    scalar E = ThermoType::Cv(Rhoest, e, T)*T; //- Initial guess

    int    iter = 0;
    do
    {
        Rhoest = Rhonew;

        E = ThermoType::Es(Rhoest, E, T);

        scalar dpdRho(-ThermoType::dpdv(Rhoest, E, T)/sqr(max(Rhoest, 1e-10)));
        Rhonew = Rhoest - (ThermoType::p(Rhoest, E, T) - p)/stabilise(dpdRho, small);

        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }
        if (Rhonew < 1e-10)
        {
            Rhonew = Rhoest/2.0;
        }
    } while (mag(Rhonew - Rhoest) > Rhotol);

    return Rhonew;
}

// ************************************************************************* //
