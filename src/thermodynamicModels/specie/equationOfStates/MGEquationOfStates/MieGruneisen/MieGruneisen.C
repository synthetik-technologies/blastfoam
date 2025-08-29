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

#include "MieGruneisen.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::MieGruneisen<Specie>::MieGruneisen(const dictionary& dict)
:
    Specie(dict),

    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    pInf_(dict.subDict("equationOfState").lookup<scalar>("pInf")),
    c0_(0.0),
    Gamma_
    (
        dict.subDict("equationOfState").found("Gamma")
      ? dict.subDict("equationOfState").lookup<scalar>("Gamma")
      : dict.subDict("equationOfState").lookup<scalar>("gamma") - 1.0
    )
{
    const dictionary& eosDict = dict.subDict("equationOfState");
    if (eosDict.found("c0"))
    {
        eosDict.readIfPresent("c0", c0_);
    }
    else if (eosDict.found("e0") && eosDict.found("p0") && eosDict.found("rho0"))
    {
        const scalar e0 = eosDict.lookup<scalar>("e0");
        const scalar p0 = eosDict.lookup<scalar>("p0");
        const scalar rho0 = eosDict.lookup<scalar>("rho0");
        c0_ = e0 - (p0 + (Gamma_ + 1.0)*pInf_)/(rho0*Gamma_);
    }
    if (Gamma_ <= 0.0)
    {
        FatalErrorInFunction
            << "gamma must be greater than 0."
            << abort(FatalError);
    }
    Info<<c0_<<endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::MieGruneisen<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("rho0", rho0_);
    dict.add("pInf", pInf_);
    dict.add("c0", c0_);
    dict.add("Gamma", Gamma_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const MieGruneisen<Specie>& ig
)
{
    ig.write(os);
    return os;
}


// ************************************************************************* //
