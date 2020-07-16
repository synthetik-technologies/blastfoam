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

#include "AUSMPlus.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(AUSMPlus, 0);
    addToRunTimeSelectionTable(fluxScheme, AUSMPlus, dictionary);
}
}


Foam::scalar Foam::fluxSchemes::AUSMPlus::f
(
    const scalar& M,
    const label sign
) const
{
    scalar magM(mag(M));

    if (magM >= 0)
    {
        return 0.5*(M + sign*(magM));
    }
    else
    {
        return sign*0.25*sqr(M + sign) + sign*0.125*sqr(sqr(M) - 1);
    }
}

Foam::scalar Foam::fluxSchemes::AUSMPlus::beta
(
    const scalar& M,
    const label sign,
    const scalar& fa
) const
{
    scalar magM(mag(M));

    if (magM >= 0)
    {
        return 0.5*(1.0 + Foam::sign(sign*M));
    }
    else
    {
        scalar A(3.0/16.0*(5.0*sqr(fa) - 4.0));
        return
            0.25*(2.0 - sign*M)*sqr(M + sign*1.0)
          + sign*A*M*sqr(sqr(M) - 1);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::AUSMPlus::AUSMPlus(const fvMesh& mesh, const word& name)
:
    fluxScheme(mesh, name),
    ku_(1.0),
    kp_(1.0),
    cutOffMa_(1e-10),
    residualRho_(1e-10),
    residualU_(1e-10)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::AUSMPlus::~AUSMPlus()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fluxSchemes::AUSMPlus::clear()
{
    fluxScheme::clear();
    phi_.clear();
}

void Foam::fluxSchemes::AUSMPlus::createSavedFields()
{
    fluxScheme::createSavedFields();
    if (phi_.valid())
    {
        return;
    }
    phi_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                IOobject::groupName("AUSMPlus::phi", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
}


void Foam::fluxSchemes::AUSMPlus::calculateFluxes
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const scalar& rhoO, const scalar& rhoN,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalar& alphaPhi,
    scalar& alphaRhoPhi,
    vector& alphaRhoUPhi,
    scalar& alphaRhoEPhi,
    const label facei, const label patchi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    scalar rhoOwn = max(rhoO, 1e-6);
    scalar rhoNei = max(rhoN, 1e-6);
    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar HOwn(EOwn + pOwn/rhoOwn);

    scalar ENei = eNei + 0.5*magSqr(UNei);
    scalar HNei(ENei + pNei/rhoNei);

    scalar UvOwn(UOwn & normal);
    scalar UvNei(UNei & normal);

    scalar cStar(sqrt(cOwn*cNei));

    // Compute split Mach numbers
    scalar MaOwn(UvOwn/cStar);
    scalar MaNei(UvNei/cStar);

    scalar MaStar(f(MaOwn, 1.0) + f(MaNei, -1.0));
    scalar M0
    (
        sqrt(min(1.0, max(sqr((MaOwn + MaNei)*0.5), cutOffMa_)))
    );
    scalar fa(M0*(2.0 - M0));

    scalar BetaOwn(beta(MaOwn, 1.0, fa));
    scalar BetaNei(beta(MaNei, -1.0, fa));

    scalar deltaMa
    (
        f(MaOwn, 1.0) - pos0(MaOwn) - f(MaNei, -1.0) + neg(MaNei)
    );
    scalar Du
    (
        -ku_*BetaOwn*BetaNei
       *0.5*(alphaOwn*rhoOwn + alphaNei*rhoNei)
       *fa*cStar*(UvNei - UvOwn)
    );
    scalar Dp
    (
        -kp_/fa*deltaMa*max(1.0 - sqr(0.5*(MaOwn - MaNei)), 0.0)
       *(alphaNei*pNei - alphaOwn*pOwn)/cStar
    );

    scalar mDot
    (
        0.5*cStar
       *(
            alphaOwn*rhoOwn*max(MaStar, 0.0)
          + alphaNei*rhoNei*min(MaStar, 0.0)
        ) + Dp
    );
    scalar alphaP
    (
        BetaOwn*alphaOwn*pOwn + BetaNei*alphaNei*pNei + Du
    );

    alphaRhoPhi = magSf*mDot;
    scalar alpha;
    vector U;
    scalar rho;

    if (mDot >= 0)
    {
        alpha = alphaOwn;
        rho = rhoOwn;
        U = save(facei, patchi, UOwn, Uf_);
    }
    else
    {
        alpha = alphaNei;
        rho = rhoNei;
        U = save(facei, patchi, UNei, Uf_);
    }
    phi = save(facei, patchi, U & Sf, phi_);
    alphaPhi = alphaRhoPhi/rho;
    save
    (
        facei,
        patchi,
        alphaP/max(alpha, 1e-6),
        pf_
    );

    alphaRhoUPhi =
        magSf*0.5
       *(
            mDot*(UOwn + UNei)
          + mag(mDot)*(UOwn - UNei)
        )
      + alphaP*Sf;

    alphaRhoEPhi =
        magSf*0.5
       *(
            mDot*(HOwn + HNei)
          + mag(mDot)*(HOwn - HNei)
        );
}


Foam::scalar Foam::fluxSchemes::AUSMPlus::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const bool rho,
    const label facei, const label patchi
) const
{

    if (getValue(facei, patchi, phi_()) >= 0)
    {
        return fOwn;
    }
    else
    {
        return fNei;
    }
}
