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

#include "VanderWaals.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::VanderWaals<Specie>::VanderWaals
(
    const dictionary& dict
)
:
    Specie(dict),
    a_(0),
    b_(0),
    Tc_(0),
    Pc_(0)
{
    const dictionary& eosDict = dict.subDict("equationOfState");
    if (eosDict.found("a") && eosDict.found("b"))
    {
        a_ = eosDict.lookup<scalar>("a");
        b_ = eosDict.lookup<scalar>("b");
        Pc_ = a_/(27.0*sqr(b_));
        Tc_ = 8*a_/(27.0*b_*this->R());
    }
    else if (eosDict.found("Tc") && eosDict.found("Pc"))
    {
        Tc_ = eosDict.lookup<scalar>("Tc");
        Pc_ = eosDict.lookup<scalar>("Pc");
        a_ = 27.0*sqr(Tc_*this->R())/(64.0*Pc_);
        b_ = this->R()*Tc_/8.0/Pc_;
    }
    else if (eosDict.found("Vc") && eosDict.found("Pc"))
    {
        scalar Vc(eosDict.lookup<scalar>("Vc"));
        Pc_ = eosDict.lookup<scalar>("Pc");
        a_ = 3.0*Pc_*sqr(Vc);
        b_ = Vc/3.0;
        Tc_ = 8*a_/(27.0*b_*this->R());
    }
    else
    {
        FatalIOErrorInFunction(eosDict)
            << "Could not determine parameters. Either a and b, Tc and Pc, or " << nl
            << "Vc and Pc must be provided" << endl
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::VanderWaals<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("a", a_);
    dict.add("b", b_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const VanderWaals<Specie>& VanderWaals
)
{
    VanderWaals.write(os);
    return os;
}


// ************************************************************************* //
