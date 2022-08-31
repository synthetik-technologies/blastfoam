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

#include "AUSMPlusUpPhaseFluxScheme.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace phaseFluxSchemes
{
    defineTypeNameAndDebug(AUSMPlusUp, 0);
    addToRunTimeSelectionTable(phaseFluxScheme, AUSMPlusUp, solid);
}
}



Foam::scalar Foam::phaseFluxSchemes::AUSMPlusUp::M1
(
    const scalar& M,
    const label sign
)
{
    return 0.5*(M + sign*mag(M));
}

Foam::scalar Foam::phaseFluxSchemes::AUSMPlusUp::M2
(
    const scalar& M,
    const label sign
)
{
    return sign*0.25*sqr(M + sign);
}

Foam::scalar Foam::phaseFluxSchemes::AUSMPlusUp::M4
(
    const scalar& M,
    const label sign,
    const scalar& beta
)
{
    scalar magM1(mag(M) - 1.0);
    if (magM1 >= 0)
    {
        return M1(M, sign);
    }
    else
    {
        return M2(M, sign)*(1.0 - sign*16.0*beta*M2(M, -sign));
    }
}

Foam::scalar Foam::phaseFluxSchemes::AUSMPlusUp::P5
(
    const scalar& M,
    const label sign,
    const scalar& xi
)
{
    scalar magM1(mag(M) - 1.0);
    if (magM1 >= 0)
    {
        return pos0(sign*M);
    }
    else
    {
        return M2(M, sign)*((sign*2.0 - M) - sign*xi*M*M2(M, -sign));
    }
}

void Foam::phaseFluxSchemes::AUSMPlusUp::postUpdate()
{
    alphapOwn_.clear();
    alphapNei_.clear();
    alphaMaxf_.clear();
    alphaMinFrictionf_.clear();
}


void Foam::phaseFluxSchemes::AUSMPlusUp::preUpdate(const volScalarField& p)
{
    const kineticTheorySystem& kt =
    (
        mesh_.lookupObject<kineticTheorySystem>
        (
            "kineticTheorySystem"
        )
    );

    const volScalarField& alphap(kt.alphap());
    word ktName;
    if (kt.polydisperse())
    {
        ktName = alphap.group();
    }
    else
    {
        ktName = this->group();
    }

    autoPtr<ReconstructionScheme<scalar>> alphapLimiter
    (
        ReconstructionScheme<scalar>::New(alphap, "alpha")
    );

    alphapOwn_ = alphapLimiter->interpolateOwn();
    alphapNei_ = alphapLimiter->interpolateNei();

    alphaMaxf_ = fvc::interpolate(kt.alphaMax());
    alphaMinFrictionf_ =
        fvc::interpolate(kt.alphaMinFriction());
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::AUSMPlusUp::AUSMPlusUp
(
    const fvMesh& mesh,
    const word& name
)
:
    phaseFluxScheme(mesh, name),
    beta_(dict_.lookupOrDefault("beta", 0.125)),
    fa_(dict_.lookupOrDefault("fa", 1.0)),
    D_(dict_.lookupOrDefault("D", 1.0)),
    cutOffMa_(dict_.lookupOrDefault("cutOffMa", 1e-10)),
    residualAlphaRho_(dict_.lookupOrDefault("residualAlphaRho", 1e-6)),
    residualC_(dict_.lookupOrDefault("residualC", 1e-3))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluxSchemes::AUSMPlusUp::~AUSMPlusUp()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::phaseFluxSchemes::AUSMPlusUp::clear()
{
    phaseFluxScheme::clear();
    phi_.clear();
}

void Foam::phaseFluxSchemes::AUSMPlusUp::createSavedFields()
{
    phaseFluxScheme::createSavedFields();
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
                IOobject::groupName("AUSMPlusUp::phi", this->group()),
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
}


void Foam::phaseFluxSchemes::AUSMPlusUp::calculateFluxes
(
    const scalar& alphaOwn, const scalar& alphaNei,
    const scalar& rhoOwn, const scalar& rhoNei,
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

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar c12
    (
        max
        (
            sqrt
            (
                max
                (
                    (alphaOwn*sqr(cOwn) + alphaNei*sqr(cNei))
                   /Foam::max(alphaOwn + alphaNei, 1e-6),
                    0.0
                )
            ),
            residualC_
        )
    );

    // Compute slpit Mach numbers
    scalar MaOwn(UvOwn/c12);
    scalar MaNei(UvNei/c12);
    scalar MaBarSqr((sqr(UvOwn) + sqr(UvNei))/(2.0*sqr(c12)));

    scalar Ma4Own(M4(MaOwn, 1, beta_));
    scalar Ma4Nei(M4(MaNei, -1, beta_));

    scalar alphaMax(getValue(facei, patchi, alphaMaxf_()));
    scalar alphaMinFriction(getValue(facei, patchi, alphaMinFrictionf_()));

    scalar alphaP =
        max
        (
            getValue(facei, patchi, alphapOwn_()),
            getValue(facei, patchi, alphapNei_())
        );

    scalar zeta
    (
        max
        (
            (alphaP - alphaMinFriction)/(alphaMax - alphaMinFriction),
            0.0
        )
    );
    scalar G(max(2.0*(1.0 - D_*sqr(zeta)), 0.0));

    scalar xi = 3.0/16.0*(5.0*sqr(fa_) - 4.0);
    scalar Kp(0.25 + 0.75*(1.0 - G/2.0));
    scalar Ku(0.75 + 0.25*(1.0 - G/2.0));
    scalar sigma(0.75*G/2.0);

    scalar Ma12
    (
        Ma4Own + Ma4Nei
      - 2.0*Kp/fa_*max(1.0 - sigma*MaBarSqr, 0.0)*(pNei - pOwn)
       /(
            max(alphaOwn*rhoOwn + alphaNei*rhoNei, residualAlphaRho_)
           *max(sqr(c12), residualC_)
        )
    );

    scalar p5Own(P5(MaOwn, 1, xi));
    scalar p5Nei(P5(MaNei, -1, xi));

    scalar pf =
        save
        (
            facei,
            patchi,
           -Ku*fa_*c12*p5Own*p5Nei
           *(alphaOwn*rhoOwn + alphaNei*rhoNei)*(UvNei - UvOwn)
          + p5Own*pOwn + p5Nei*pNei,
            pf_
        );
    scalar F
    (
        c12
       *(1.0 + mag(Ma12)*(1.0 - G/2.0))
       *alphaP
       *(alphaOwn*rhoOwn - alphaNei*rhoNei)/(2.0*alphaMax)
    );

    alphaRhoPhi =
        (
            F
          + c12*Ma12
           *(
                pos(Ma12)*alphaOwn*rhoOwn
              + neg0(Ma12)*alphaNei*rhoNei
            )
        )*magSf;

    scalar alpha;
    scalar rho;
    vector U;
    scalar e;
    if (alphaRhoPhi >= 0)
    {
        alpha = alphaOwn;
        rho = rhoOwn;
        U = UOwn;
        e = eOwn;
    }
    else
    {
        alpha = alphaNei;
        rho = rhoNei;
        U = UNei;
        e = eNei;
    }

    this->save(facei, patchi, alpha, alphaf_);
    this->save(facei, patchi, U, Uf_);
    this->save(facei, patchi, pf, pf_);

    phi = save
    (
        facei,
        patchi,
        alphaRhoPhi/max(alpha, 1e-6)/rho,
        phi_
    );
    alphaPhi = alphaRhoPhi/rho;
    alphaRhoUPhi = alphaRhoPhi*U + pf*Sf;
    alphaRhoEPhi = alphaRhoPhi*e;
}


Foam::scalar Foam::phaseFluxSchemes::AUSMPlusUp::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const bool rho,
    const label facei, const label patchi
) const
{

    if (getValue(facei, patchi, phi_) >= 0)
    {
        return fOwn;
    }
    else
    {
        return fNei;
    }
}
