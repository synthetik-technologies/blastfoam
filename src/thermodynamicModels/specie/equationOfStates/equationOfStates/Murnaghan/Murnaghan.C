/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

\*---------------------------------------------------------------------------*/

#include "Murnaghan.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::Murnaghan<Specie>::Murnaghan
(
    const dictionary& dict
)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookupType<scalar>("rho0")),
    pRef_(dict.subDict("equationOfState").lookupType<scalar>("pRef")),
    n_(0.0),
    kappa_(0.0),
    Gamma_(dict.subDict("equationOfState").lookupType<scalar>("Gamma"))
{
    const dictionary& eosDict = dict.subDict("equationOfState");
    if (eosDict.found("kappa"))
    {
        kappa_ = eosDict.lookupType<scalar>("kappa");
    }
    else
    {
        kappa_ = 1.0/max(eosDict.lookupType<scalar>("K0"), small);
    }

    if (eosDict.found("n"))
    {
        n_ = eosDict.lookupType<scalar>("n");
    }
    else
    {
        n_ = eosDict.lookupOrDefault<scalar>("K0Prime", 1.0);
    }
}

// ************************************************************************* //
