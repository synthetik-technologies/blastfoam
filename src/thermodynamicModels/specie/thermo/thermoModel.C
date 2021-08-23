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
    tolerance_(dict.lookupOrDefault("tolerance", 1e-4)),
    maxIter_(dict.lookupOrDefault("maxIter", 100))
{
    this->set(*this);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::initializeEnergy
(
    const scalar p,
    const scalar rho,
    const scalar e,
    const scalar T
) const
{
    if (rho < small)
    {
        return 0.0;
    }

    scalar Eest = ThermoType::Es(rho, e, T);

    if (ThermoType::temperatureBased())
    {
        return Eest;
    }

    if (mag(ThermoType::dpde(rho, e, T)) < small)
    {
        return Eest;
    }

    scalar Enew = Eest;
    scalar pNew;
    int    iter = 0;
    do
    {
        Eest = Enew;
        pNew = ThermoType::p(rho, Eest, T, false); // Do not limit p
        Enew -=
            (pNew - p)/stabilise(ThermoType::dpde(rho, Eest, T), small);

        if (iter++ > 100)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << 100
                << abort(FatalError);
        }

    } while ((mag(pNew - p)/p > tolerance_) && (iter++ < maxIter_));
    return Enew;
}


template<class ThermoType>
Foam::scalar Foam::thermoModel<ThermoType>::rhoPT
(
    const scalar p,
    const scalar T
) const
{
    //- Simple method to calculate initial density
    //  Should be modified to solve the 2D problem for
    //  density and internal energy
    scalar Rhoest = 1.0;
    scalar Rhonew = Rhoest;
    scalar pNew = p;
    scalar E = ThermoType::Es(Rhoest, 0.0, T); //- Initial guess

    int    iter = 0;
    do
    {
        Rhoest = Rhonew;

        scalar dpdRho(-ThermoType::dpdv(Rhoest, E, T)/sqr(max(Rhoest, 1e-10)));
        pNew = ThermoType::p(Rhoest, E, T, false); // Do not limit p
        Rhonew =
            Rhoest - (pNew - p)/stabilise(dpdRho, small);
        E = ThermoType::Es(Rhoest, E, T);

        if (Rhonew < 1e-10)
        {
            Rhonew = Rhoest/2.0;
        }
    } while ((mag(pNew - p)/p > tolerance_) && (iter++ < maxIter_));

    return Rhonew;
}


// ************************************************************************* //
