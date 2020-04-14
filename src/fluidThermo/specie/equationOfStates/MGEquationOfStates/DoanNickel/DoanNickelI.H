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

template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    return -dGammadRho(rho, e)*e*rho/(Gamma(rho, e) - 1.0);
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::Gamma
(
    const scalar& rho,
    const scalar& eIn
) const
{
    //- Scale internal energy to 10^10 ergs/gram
    scalar e = eIn*1e-6;
    scalar rhoByRho0(rho/rho0_);

    scalar I1(E1Offset_ + E11_*log10(rhoByRho0));
    scalar deltaI1(max(deltaE1_*pow(rhoByRho0, deltaE1Pow_), 1e-10));
    scalar f1((e - I1)/deltaI1);
    f1 = 1.0/(exp(min(f1, expMax_)) + 1.0);


    scalar I2(E22_*pow(rhoByRho0, E2Pow_));
    scalar deltaI2(max(deltaE2_*pow(rhoByRho0, deltaE2Pow_), 1e-10));
    scalar f2((I2 - e)/deltaI2);
    f2 = 1.0/(exp(min(f2, expMax_)) + 1.0);

    scalar f3((E33_ - e)/deltaE3_);
    f3 = 1.0/(exp(min(f3, expMax_)) + 1.0);

    scalar A
    (
        (a1_*f1 + a2_*(1.0 - f1)*(1.0 - f2))*log10(max(e, small))
      + a3_*f2
    );

    return
        (
            a_
          + b_*exp(-e/E1_)*f1
          + c_*exp(-e/E2_)*(1.0 - f1)
          + d_*exp(-e/E3_)*f2
          + g_*f3
        )*pow(rhoByRho0, A) + 1.0;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    return -sqr(rho)*(dGammadRho(rho, e)*rho*e + (Gamma(rho, e) - 1.0)*e);
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    //- Unscale the change in internal energy
    return 1e-6*dGammade(rho, e)*rho*e + (Gamma(rho, e) - 1.0)*rho;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    // Need to calculate
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    // Need to calculate
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    // Need to calculate
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    // Need to calculate
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::CpMCv
(
    const scalar& rho,
    const scalar&e
) const
{
    // Need to calculate
    return this->R();
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    // Need to calculate
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::dGammadRho
(
    const scalar& rho,
    const scalar& eIn
) const
{
    scalar e = eIn*1e-6;
    scalar rhoByRho0(rho/rho0_);

    scalar I1(E1Offset_ + E11_*log10(rhoByRho0));
    scalar deltaI1(max(deltaE1_*pow(rhoByRho0, deltaE1Pow_), 1e-10));
    scalar f1((e - I1)/deltaI1);
    scalar df1 = 0.0;
    if (f1 > expMax_)
    {
        f1 = 0.0;
    }
    else
    {
        scalar dI1(E11_/(max(rho, 1e-10)*log(10.0)));
        scalar dDeltaI1(deltaE2Pow_*deltaE2_*pow(rhoByRho0, deltaE2Pow_ - 1.0));

        df1 =
            exp(f1)/sqr(exp(f1) + 1.0)
           *(dI1/deltaI1 + dDeltaI1*(e - I1)/sqr(deltaI1));
        f1 = 1.0/(exp(f1) + 1.0);
    }

    scalar I2(E22_*pow(rhoByRho0, E2Pow_));
    scalar deltaI2(max(deltaE2_*pow(rhoByRho0, deltaE2Pow_), 1e-10));
    scalar f2((I2 - e)/deltaI2);
    scalar df2 = 0.0;
    if (f2 > expMax_)
    {
        f2 = 0.0;
    }
    else
    {
        scalar dI2(E22_*E2Pow_/rho0_*pow(rhoByRho0, E2Pow_ - 1.0));
        scalar dDeltaI2(deltaE2Pow_*deltaE2_/rho0_*pow(rhoByRho0, deltaE2Pow_ - 1.0));

        df2 =
          - exp(f2)/sqr(exp(f2) + 1.0)
           *(dI2/deltaI2 - dDeltaI2*(I2 - e)/sqr(deltaI2));
        f2 = 1.0/(exp(f2) + 1.0);
    }

    scalar f3((E33_ - e)/deltaE3_);
    if (f3 > expMax_)
    {
        f3 = 0.0;
    }
    else
    {
        f3 = 1.0/(exp(f3) + 1.0);
    }

    scalar A
    (
        (a1_*f1 + a2_*(1.0 - f1)*(1.0 - f2))*log10(max(e, small))
      + a3_*f2
    );

    scalar dA
    (
        (a1_*df1 - a2_*(df1*(1.0 - f2) + (1.0 - f1)*df2))*log10(max(e, small))
      + a3_*df2
    );

    scalar g =
            a_
          + b_*exp(-e/E1_)*f1
          + c_*exp(-e/E2_)*(1.0 - f1)
          + d_*exp(-e/E3_)*f2
          + g_*f3;
    scalar dg =
            df1*b_*exp(-e/E1_)
          - df1*c_*exp(-e/E2_)
          + df2*d_*exp(-e/E3_);

    return (g*(log(rhoByRho0)*dA + A/rho) + dg)*pow(rhoByRho0, A);
}


template<class Specie>
Foam::scalar Foam::DoanNickel<Specie>::dGammade
(
    const scalar& rho,
    const scalar& eIn
) const
{
    scalar e = eIn*1e-6;
    scalar rhoByRho0(rho/rho0_);

    scalar I1(E1Offset_ + E11_*log10(rhoByRho0));
    scalar deltaI1(max(deltaE1_*pow(rhoByRho0, deltaE1Pow_), 1e-10));
    scalar f1((e - I1)/deltaI1);
    scalar df1 = 0.0;
    if (f1 > expMax_)
    {
        f1 = 0.0;
    }
    else
    {
        df1 = -exp(f1)/(sqr(exp(f1) + 1.0)*deltaI1);
        f1 = 1.0/(exp(f1) + 1.0);
    }

    scalar I2(E22_*pow(rhoByRho0, E2Pow_));
    scalar deltaI2(max(deltaE2_*pow(rhoByRho0, deltaE2Pow_), 1e-10));
    scalar f2((I2 - e)/deltaI2);
    scalar df2 = 0.0;
    if (f2 > expMax_)
    {
        f2 = 0.0;
    }
    else
    {
        df2 = exp(f2)/(sqr(exp(f2) + 1.0)*deltaI2);
        f2 = 1.0/(exp(f2) + 1.0);
    }

    scalar f3((E33_ - e)/deltaE3_);
    scalar df3 = 0.0;
    if (f3 > expMax_)
    {
        f3 = 0.0;
    }
    else
    {
        df3 = exp(f3)/(sqr(exp(f3) + 1.0)*deltaE3_);
        f3 = 1.0/(exp(f3) + 1.0);
    }

    scalar A
    (
        (a1_*f1 + a2_*(1.0 - f1)*(1.0 - f2))*log10(max(e, small))
      + a3_*f2
    );

    scalar dA
    (
        (a1_*df1 - a2_*(df1*(1.0 - f2) + (1.0 - f1)*df2))*log10(max(e, small))
      + (a1_*f1 + a2_*(1.0 - f1)*(1.0 - f2))/(e*log(10.0))
      + a3_*df2
    );

    scalar g =
            a_
          + b_*exp(-e/E1_)*f1
          + c_*exp(-e/E2_)*(1.0 - f1)
          + d_*exp(-e/E3_)*f2
          + g_*f3;
    scalar dg =
            b_*exp(-e/E1_)*(df1 - f1/E1_)
          - c_*exp(-e/E2_)*(df1 + (1.0 - f1)/E2_)
          + d_*exp(-e/E3_)*(df2 - f2/E3_)
          + g_*df3;

    return (g*log(rhoByRho0)*dA + dg)*pow(rhoByRho0, A);
}

// ************************************************************************* //
