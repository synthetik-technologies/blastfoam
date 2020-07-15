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

#include "Tillotson.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::Tillotson<Specie>::Tillotson(const dictionary& dict)
:
    Specie(dict),
    p0_(dict.subDict("equationOfState").lookupType<scalar>("p0")),
    pCav_(dict.subDict("equationOfState").lookupType<scalar>("pCav")),
    rho0_(dict.subDict("equationOfState").lookupType<scalar>("rho0")),
    rhoCav_(dict.subDict("equationOfState").lookupType<scalar>("rhoCav")),
    e0_(dict.subDict("equationOfState").lookupType<scalar>("e0")),
    omega_(dict.subDict("equationOfState").lookupType<scalar>("omega")),
    k_(dict.subDict("equationOfState").lookupType<scalar>("k")),
    A_(dict.subDict("equationOfState").lookupType<scalar>("A")),
    B_(dict.subDict("equationOfState").lookupType<scalar>("B")),
    C_(dict.subDict("equationOfState").lookupType<scalar>("C")),
    EcTable_()
{
    //- Dummy temperature
    scalar T(273.0);

    label tableSize(100);
    scalarField rhof(tableSize + 1, rho0_);
    scalarField ecf(tableSize + 1, 0.0);
    label I = 1;

    scalar x = rho0_;
    scalar y = 0.0;
    label nSteps = 10000;
    scalar dx = 0.5*rho0_/scalar(nSteps);


    for (label i = 0; i < nSteps; i++)
    {
        scalar yOld = y;

        scalar k1 =
            (((Gamma(x, y, T, T) - 1.0)*x*y) - Pi(x, y))/sqr(x);
        y = yOld + dx*0.5*k1;

        x += 0.5*dx;
        scalar k2 =
            (((Gamma(x, y, T, T) - 1.0)*x*y) - Pi(x, y))/sqr(x);
        y = yOld + dx*0.5*k2;
        scalar k3 =
            (((Gamma(x, y, T, T) - 1.0)*x*y) - Pi(x, y))/sqr(x);

        y = yOld + dx*k3;
        x += 0.5*dx;
        scalar k4 =
            (((Gamma(x, y, T, T) - 1.0)*x*y) - Pi(x, y))/sqr(x);

        y = yOld + dx/6.0*(k1 + 2.0*(k2 + k3) + k4);

        if (((i + 1) % (nSteps/tableSize)) == 0)
        {
            rhof[I] = x;
            ecf[I] = y;
            I++;
        }
    }

    EcTable_.set
    (
        rhof,
        ecf,
        "none",
        "none",
        true
    );
}

// ************************************************************************* //
