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

#include "PengRobinson.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::PengRobinson<Specie>::PengRobinson
(
    const dictionary& dict
)
:
    Specie(dict),
    Tc_(dict.subDict("equationOfState").lookup<scalar>("Tc")),
    Vc_(dict.subDict("equationOfState").lookup<scalar>("Vc")),
    Pc_(dict.subDict("equationOfState").lookup<scalar>("Pc")),
    omega_(dict.subDict("equationOfState").lookup<scalar>("omega")),
    Zc_(Pc_*Vc_/(RR*Tc_)),
    a_(0.45724*sqr(RR*Tc_)/Pc_),
    b_(0.07780*RR*Tc_/Pc_),
    kappa_(0.37464 + 1.54226*omega_ - 0.26992*sqr(omega_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::PengRobinson<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("Tc", Tc_);
    dict.add("Vc", Vc_);
    dict.add("Pc", Pc_);
    dict.add("omega", omega_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const PengRobinson<Specie>& pr
)
{
    pr.write(os);
    return os;
}


// ************************************************************************* //
