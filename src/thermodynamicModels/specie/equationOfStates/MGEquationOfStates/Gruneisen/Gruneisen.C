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

#include "Gruneisen.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::Gruneisen<Specie>::Gruneisen(const dictionary& dict)
:
    Specie(dict),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    C_(dict.subDict("equationOfState").lookup<scalar>("C")),
    S1_(dict.subDict("equationOfState").lookup<scalar>("S1")),
    S2_(dict.subDict("equationOfState").lookup<scalar>("S2")),
    S3_(dict.subDict("equationOfState").lookup<scalar>("S3")),
    Gamma_(dict.subDict("equationOfState").lookup<scalar>("Gamma")),
    a_(dict.subDict("equationOfState").lookup<scalar>("a"))
{
    if (Gamma_ <= 0.0)
    {
        FatalErrorInFunction
            << "gamma must be greater than 0."
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::Gruneisen<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("rho0", rho0_);
    dict.add("C", C_);
    dict.add("S1", S1_);
    dict.add("S2", S2_);
    dict.add("S3", S3_);
    dict.add("Gamma", Gamma_);
    dict.add("a", a_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Gruneisen<Specie>& ig
)
{
    ig.write(os);
    return os;
}


// ************************************************************************* //
