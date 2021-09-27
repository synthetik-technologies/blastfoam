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

#include "linearTillotson.H"
#include "lookupTable1D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::linearTillotson<Specie>::linearTillotson(const dictionary& dict)
:
    Specie(dict),
    p0_(dict.subDict("equationOfState").lookup<scalar>("p0")),
    rho0_(dict.subDict("equationOfState").lookup<scalar>("rho0")),
    e0_(dict.subDict("equationOfState").lookup<scalar>("e0")),
    omega_(dict.subDict("equationOfState").lookup<scalar>("omega")),
    A_(dict.subDict("equationOfState").lookup<scalar>("A")),
    B_(dict.subDict("equationOfState").lookup<scalar>("B")),
    C_(dict.subDict("equationOfState").lookup<scalar>("C")),
    pCav_(dict.subDict("equationOfState").lookup<scalar>("pCav"))
{
    scalar Cv = dict.subDict("thermodynamics").lookup<scalar>("Cv");
    scalar T = e0_/Cv;

    label tableSize(100);
    scalarField rhof(tableSize, rho0_);
    scalarField ecf(tableSize, 0.0);

    scalar rhoMin = dict.subDict("equationOfState").lookupOrDefault<scalar>
    (
        "rhoMin",
        0.5*rho0_
    );
    scalar rhoMax = dict.subDict("equationOfState").lookupOrDefault<scalar>
    (
        "rhoMax",
        5.0*rho0_
    );

    label n = 10000;
    scalar dRho = (rhoMax - rhoMin)/scalar(n);
    label n1 = (rho0_ - rhoMin)/dRho;
    label n2 = n - n1;

    label stepsPerEntry(n/tableSize);

    label I = n1/stepsPerEntry+1;
    scalar rho = rho0_;
    scalar ec = 0.0;
    for (label i = 0; i < n2; i++)
    {
        scalar ecOld = ec;

        scalar k1 = p(rho, ec, T)/sqr(rho);
        ec = ecOld + dRho*0.5*k1;

        rho += 0.5*dRho;
        scalar k2 = p(rho, ec, T)/sqr(rho);
        ec = ecOld + dRho*0.5*k2;
        scalar k3 = p(rho, ec, T)/sqr(rho);

        ec = ecOld + dRho*k3;
        rho += 0.5*dRho;
        scalar k4 = p(rho, ec, T)/sqr(rho);

        ec = ecOld + dRho/6.0*(k1 + 2.0*(k2 + k3) + k4);

        if (((i + 1) % label(n/tableSize)) == 0)
        {
            rhof[I] = rho;
            ecf[I] = ec;
            I++;
        }
    }

    I = n1/stepsPerEntry;
    rho = rho0_;
    ec = 0.0;
    dRho *= -1;
    for (label i = 0; i < n1; i++)
    {
        scalar ecOld = ec;

        scalar k1 = p(rho, ec, T)/sqr(rho);
        ec = ecOld + dRho*0.5*k1;

        rho += 0.5*dRho;
        scalar k2 = p(rho, ec, T)/sqr(rho);
        ec = ecOld + dRho*0.5*k2;
        scalar k3 = p(rho, ec, T)/sqr(rho);

        ec = ecOld + dRho*k3;
        rho += 0.5*dRho;
        scalar k4 = p(rho, ec, T)/sqr(rho);

        ec = ecOld + dRho/6.0*(k1 + 2.0*(k2 + k3) + k4);

        if (((i - 1) % (n/tableSize)) == 0)
        {
            rhof[I] = rho;
            ecf[I] = ec;
            I--;
        }
    }

    EcTable_.set
    (
        rhof,
        ecf,
        "none",
        "none",
        "linearExtrapolated",
        false
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class Specie>
void Foam::linearTillotson<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("p0", p0_);
    dict.add("rho0", rho0_);
    dict.add("omega", omega_);
    dict.add("e0", e0_);
    dict.add("A", A_);
    dict.add("B", B_);
    dict.add("C", C_);
    dict.add("pCav", pCav_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const linearTillotson<Specie>& lt
)
{
    lt.write(os);
    return os;
}


// ************************************************************************* //
