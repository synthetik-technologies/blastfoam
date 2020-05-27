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
Foam::scalar Foam::Tillotson<Specie>::Pi1
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return -(A_*mu + B_*sqr(mu));
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Pi2
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return
      - exp(-alpha_*sqr(rho0_/max(rho, 1e-10) - 1.0))
       *(
            b_*rho*e/(e/(e0_*sqr(eta)) + 1.0)
          + A_*mu*exp(-beta_*(rho0_/max(rho, 1e-10) - 1.0))
        );
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Pi3
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return -A_*mu;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dPidRho1
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return -(A_ + 2.0*B_*mu)/rho0_;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dPidRho2
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    scalar rhos(max(rho, 1e-10));
    scalar A = b_*rho*e/(e/(e0_*sqr(eta)) + 1.0);
    scalar APrime =
        b_*e*e0_*sqr(rho)*(3.0*e*sqr(rho0_) + e0_*sqr(rho))
       /(e*sqr(rho0_) + e0_*sqr(rho));

    scalar B = A_*mu*exp(-beta_*(rho0_/max(rho, 1e-10) - 1.0));
    scalar BPrime =
        A_*(beta_*rho0_*(rho - rho0_) + sqr(rho))
        *exp(beta_*(rho - rho0_)/sqr(rhos))
       /(sqr(rhos*rho0_));

    scalar C = exp(-alpha_*sqr(rho0_/max(rho, 1e-10) - 1.0));
    scalar CPrime =
        2.0*alpha_*rho0_*(rho - rho0_)
        *exp(-alpha_*sqr((rho - rho0_)/rhos))
       /pow3(rhos);

    return CPrime*(A + B) + C*(APrime + BPrime);
}

template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dPidRho3
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return -A_/rho0_;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dPide2
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return
        -b_*sqr(e0_)*pow5(rho)*exp(-alpha_*sqr(rho - rho0_)/sqr(max(rho, 1e-10)))
       /sqr(e*sqr(rho0_) + e0_*sqr(rho));
}

template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Gamma1
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return a_ + b_/(e/(e0_*sqr(eta)) + 1.0) + 1.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Gamma2
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return a_ + 1.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dGammadRho1
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    scalar rhos(max(rho, 1e-10));
    return
        2.0*b_*e*sqr(rho)
       /(e0_*pow3(rhos)*sqr(e*sqr(rho0_)/(e0_*sqr(rhos)) + 1.0));
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dGammade1
(
    const scalar& rho,
    const scalar& e,
    const scalar& eta,
    const scalar& mu
) const
{
    return -b_*sqr(rho0_*rho)*e0_/sqr(e*sqr(rho0_) + e0_*sqr(rho));
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Pi
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar eta(rho/rho0_);
    scalar mu(eta - 1.0);

    if
    (
        (rho >= rho0_ && e >= 0)
     || (rho >= rhoIV_ &&  rho < rho0_ && e <= eIV_)
    )
    {
        return Pi1(rho, e, eta, mu);
    }
    else if (rho < rho0_ && e >= eCV_)
    {
        return Pi2(rho, e, eta, mu);
    }
    else if (rho < rho0_ && rho > rhoIV_ && e > eIV_ && e < eCV_)
    {
        return
            ((e - eIV_)*Pi2(rho, e, eta, mu) + (eCV_ - e)*Pi1(rho, e, eta, mu))
           /(eCV_ - eIV_);
    }

    return Pi3(rho, e, eta, mu);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::delta
(
    const scalar& p,
    const scalar& rho,
    const scalar& e
) const
{
    scalar eta(rho/rho0_);
    scalar mu(eta - 1.0);
    scalar pi, dPi, gamma, dGamma;
    if
    (
        (rho >= rho0_ && e >= 0)
     || (rho >= rhoIV_ &&  rho < rho0_ && e <= eIV_)
    )
    {
        pi = Pi1(rho, e, eta, mu);
        dPi = dPidRho1(rho, e, eta, mu);
        gamma = Gamma1(rho, e, eta, mu);
        dGamma = dGammadRho1(rho, e, eta, mu);
    }
    else if ((rho < rho0_ && e >= eCV_))
    {
        pi = Pi2(rho, e, eta, mu);
        dPi = dPidRho2(rho, e, eta, mu);
        gamma = Gamma2(rho, e, eta, mu);
        dGamma = 0.0;
    }
    else if (rho < rho0_ && rho > rhoIV_ && e > eIV_ && e < eCV_)
    {
        scalar f((e - eIV_)/(eCV_ - eIV_));
        pi = f*Pi2(rho, e, eta, mu) + (1.0 - f)*Pi1(rho, e, eta, mu);
        dPi = f*dPidRho2(rho, e, eta, mu) + (1.0 - f)*dPidRho1(rho, e, eta, mu);
        gamma = f*Gamma2(rho, e, eta, mu) + (1.0 - f)*Gamma1(rho, e, eta, mu);
        dGamma = (1.0 - f)*dGammadRho1(rho, e, eta, mu);
    }
    else
    {
        pi = Pi3(rho, e, eta, mu);
        dPi = dPidRho3(rho, e, eta, mu);
        gamma = Gamma1(rho, e, eta, mu);
        dGamma = dGammadRho1(rho, e, eta, mu);
    }
    return -(p + pi)*dGamma/sqr(gamma - 1.0) + dPi/(gamma - 1.0);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Gamma
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar eta(rho/rho0_);
    scalar mu(eta - 1.0);

    if
    (
        (rho >= rhoIV_ && e <= eIV_)
     || (rho < rhoIV_ && e < eCV_)
    )
    {
        return Gamma1(rho, e, eta, mu);
    }
    else if (rho < rho0_ && e >= eCV_)
    {
        return Gamma2(rho, e, eta, mu);
    }
    return
        ((e - eIV_)*Gamma2(rho, e, eta, mu) + (eCV_ - e)*Gamma1(rho, e, eta, mu))
       /(eCV_ - eIV_);
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpdv
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar eta(rho/rho0_);
    scalar mu(eta - 1.0);
    scalar dPi, gamma, dGamma;
    if (rho >= rhoIV_ && e <= eIV_)
    {
        dPi = dPidRho1(rho, e, eta, mu);
        gamma = Gamma1(rho, e, eta, mu);
        dGamma = dGammadRho1(rho, e, eta, mu);
    }
    else if ((rho < rho0_ && e >= eCV_))
    {
        dPi = dPidRho2(rho, e, eta, mu);
        gamma = Gamma2(rho, e, eta, mu);
        dGamma = 0.0;
    }
    else if (rho < rho0_ && rho > rhoIV_ && e > eIV_ && e < eCV_)
    {
        scalar f((e - eIV_)/(eCV_ - eIV_));
        dPi = f*dPidRho2(rho, e, eta, mu) + (1.0 - f)*dPidRho1(rho, e, eta, mu);
        gamma = f*Gamma2(rho, e, eta, mu) + (1.0 - f)*Gamma1(rho, e, eta, mu);
        dGamma = (1.0 - f)*dGammadRho1(rho, e, eta, mu);
    }
    else
    {
        dPi = dPidRho3(rho, e, eta, mu);
        gamma = Gamma1(rho, e, eta, mu);
        dGamma = dGammadRho1(rho, e, eta, mu);
    }
    scalar dpdrho = dGamma*rho*e - (gamma - 1.0)*e - dPi;
    return -sqr(rho)*dpdrho;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::dpde
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar eta(rho/rho0_);
    scalar mu(eta - 1.0);
    scalar dPi, gamma, dGamma;
    if (rho >= rhoIV_ && e <= eIV_)
    {
        dPi = 0.0;
        gamma = Gamma1(rho, e, eta, mu);
        dGamma = dGammade1(rho, e, eta, mu);
    }
    else if ((rho < rho0_ && e >= eCV_))
    {
        dPi = dPide2(rho, e, eta, mu);
        gamma = Gamma2(rho, e, eta, mu);
        dGamma = 0.0;
    }
    else if (rho < rho0_ && rho > rhoIV_ && e > eIV_ && e < eCV_)
    {
        scalar f((e - eIV_)/(eCV_ - eIV_));
        dPi = f*dPide2(rho, e, eta, mu);
        gamma = f*Gamma2(rho, e, eta, mu) + (1.0 - f)*Gamma1(rho, e, eta, mu);
        dGamma = (1.0 - f)*dGammade1(rho, e, eta, mu);
    }
    else
    {
        dPi = 0.0;
        gamma = Gamma1(rho, e, eta, mu);
        dGamma = dGammade1(rho, e, eta, mu);
    }
    return dGamma*rho*e + (gamma - 1.0)*rho - dPi;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::E
(
    const scalar& rho,
    const scalar& e
) const
{
    scalar rhoOld = rho0_;

    label nSteps = 10;
    scalar Ec = 0.0;
    scalar Ec0 = 0.0;
    scalar f1, f2, f3, f4;
    scalar dRho = (rho - rhoOld)/scalar(nSteps + 1);
    scalar rho0 = rhoOld;
    scalar rho1 = rhoOld;
    scalar rho2 = rhoOld;
    scalar rho3 = rhoOld;
    scalar Ec1, Ec2, Ec3;
    for (label stepi = 0; stepi <= nSteps; stepi++)
    {
        Ec0 = Ec;
        rho0 = rho3;
        rho1 = rho0 + dRho/2.0;
        rho2 = rho1;
        rho3 = rho0 + dRho;

        f1 =
            ((Gamma(rho0, Ec0) - 1.0)*rho0*Ec0 - Pi(rho0, Ec0))/sqr(max(rho0, 1e-10));
        Ec1 = Ec0 + f1*dRho/2.0;

        f2 =
            ((Gamma(rho1, Ec1) - 1.0)*rho1*Ec1 - Pi(rho1, Ec1))/sqr(max(rho1, 1e-10));
        Ec2 = Ec1 + f2*dRho/2.0;

        f3 =
            ((Gamma(rho2, Ec2) - 1.0)*rho2*Ec2 - Pi(rho2, Ec2))/sqr(max(rho2, 1e-10));
        Ec3 = Ec0+ f3*dRho;

        f4 =
            ((Gamma(rho3, Ec3) - 1.0)*rho3*Ec3 - Pi(rho3, Ec3))/sqr(max(rho3, 1e-10));

        Ec = Ec0 + 1.0/6.0*dRho*(f1 + 2.0*(f2 + f3) + f4);
    }
    return Ec;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Cv
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::Cp
(
    const scalar& rho,
    const scalar& e
) const
{
    // Check
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::H
(
    const scalar& rho,
    const scalar& e
) const
{
    return 0.0;
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::CpMCv
(
    const scalar& rho,
    const scalar& e
) const
{
    return this->R();
}


template<class Specie>
Foam::scalar Foam::Tillotson<Specie>::S
(
    const scalar& p,
    const scalar& rho,
    const scalar& T
) const
{
    return -this->R()*log(p/Foam::constant::thermodynamic::Pstd);
}

// ************************************************************************* //
