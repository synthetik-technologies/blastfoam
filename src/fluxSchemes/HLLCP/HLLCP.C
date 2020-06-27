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

#include "HLLCP.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fluxSchemes
{
    defineTypeNameAndDebug(HLLCP, 0);
    addToRunTimeSelectionTable(fluxScheme, HLLCP, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLLCP::HLLCP
(
    const fvMesh& mesh
)
:
    fluxScheme(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluxSchemes::HLLCP::~HLLCP()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluxSchemes::HLLCP::clear()
{
    fluxScheme::clear();
    SOwn_.clear();
    SNei_.clear();
    SStar_.clear();
    UvOwn_.clear();
    UvNei_.clear();
}

void Foam::fluxSchemes::HLLCP::createSavedFields()
{
    fluxScheme::createSavedFields();
    if (SOwn_.valid())
    {
        return;
    }
    SOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::SOwn",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    SNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::SNei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    SStar_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::SStar",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    pStar_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::pStar",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimPressure, 0.0)
        )
    );
    phip_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::phip",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity*dimArea, 0.0)
        )
    );
    UTilde_ = tmp<surfaceVectorField>
    (
        new surfaceVectorField
        (
            IOobject
            (
                "HLLCP::UTilde",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedVector("0", dimVelocity, Zero)
        )
    );
    UvOwn_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::UvOwn",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
    UvNei_ = tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                "HLLCP::UvNei",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("0", dimVelocity, 0.0)
        )
    );
}

void Foam::fluxSchemes::HLLCP::preUpdate(const volScalarField& p)
{
    volScalarField fCells
    (
        IOobject
        (
            "fCells",
            p.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar("0", dimless, 1.0)
    );

    const labelListList& cellCells = mesh_.cellCells();
    forAll(fCells, celli)
    {
        forAll(cellCells[celli], cellj)
        {
            fCells[celli] =
                min
                (
                    fCells[celli],
                    min
                    (
                        p[celli]/p[cellCells[celli][cellj]],
                        p[cellCells[celli][cellj]]/p[celli]
                    )
                );
        }
    }
    f_ = fvc::interpolate(pow3(fCells));
}


void Foam::fluxSchemes::HLLCP::postUpdate()
{
    f_.clear();
}


void Foam::fluxSchemes::HLLCP::calculateFluxes
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalar& rhoPhi,
    vector& rhoUPhi,
    scalar& rhoEPhi,
    const label facei, const label patchi
)
{
    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar MaOwn(UvOwn/cOwn);
    scalar MaNei(UvNei/cNei);

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(sqrt(rhoNei)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar MaTilde(wOwn*MaOwn + wNei*MaNei);
    vector UTilde(wOwn*UOwn + wNei*UNei);
    scalar UvTilde(wOwn*UvOwn + wNei*UvNei);
    scalar cTilde(wOwn*cOwn + wNei*cNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    scalar aOwn(rhoOwn*(SOwn - UvOwn));
    scalar aNei(rhoNei*(SNei - UvNei));

    scalar SStar((aNei*UvNei - aOwn*UvOwn + pOwn - pNei)/(aNei - aOwn));
    scalar pStar
    (
        (
            aNei*pOwn - aOwn*pNei - aOwn*aNei*(UvOwn - UvNei)
        )/(aNei - aOwn)
    );

    scalar theta(min(max(mag(MaOwn), mag(MaNei)), 1.0));
    scalar pAvg(0.5*(pOwn + pNei));
    scalar pStarStar(pStar*theta + (1.0 - theta)*pAvg);

    scalar f(f_()[facei]);
    scalar pStarStarStar(f*pStarStar + (1.0 - f)*pStar);

    scalar phip =
        (f - 1.0)
       *SOwn*SNei/(SNei - SOwn)
       /(1.0 + mag(MaTilde))
       *(pNei - pOwn)/sqr(cTilde);

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);
    this->save(facei, patchi, pStar, pStar_);
    this->save(facei, patchi, phip, phip_);
    this->save(facei, patchi, UTilde, UTilde_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    const vector rhoUOwn = rhoOwn*UOwn;
    const scalar rhoEOwn = rhoOwn*EOwn;

    const vector rhoUPhiOwn = rhoUOwn*UvOwn + pOwn*normal;
    const scalar rhoEPhiOwn = (rhoEOwn + pOwn)*UvOwn;

    // Neighbour values
    const vector rhoUNei = rhoNei*UNei;
    const scalar rhoENei = rhoNei*ENei;

    const vector rhoUPhiNei = rhoUNei*UvNei + pNei*normal;
    const scalar rhoEPhiNei = (rhoENei + pNei)*UvNei;

    scalar p;
    if (SOwn > 0)
    {
        this->save(facei, patchi, UOwn, Uf_);
        phi = UvOwn;
        rhoPhi = rhoOwn*UvOwn;
        rhoUPhi = rhoUPhiOwn;
        rhoEPhi = rhoEPhiOwn;
        p = pOwn;
    }
    else if (SStar > 0)
    {
        const scalar dS = SOwn - SStar;

        this->save
        (
            facei,
            patchi,
            (SOwn*rhoUOwn - rhoUPhiOwn + pStarStarStar*normal)
           /(rhoOwn*(SOwn - UvOwn)),
            Uf_
        );
        phi = SStar*(SOwn - UvOwn)/dS;

        rhoPhi = phi*rhoOwn + phip;
        rhoUPhi =
            (
                SStar*(SOwn*rhoUOwn - rhoUPhiOwn) + SOwn*pStarStarStar*normal
            )/dS
          + phip*UTilde;
        rhoEPhi =
            SStar*(SOwn*rhoEOwn - rhoEPhiOwn + SOwn*pStar)/dS
          + 0.5*phip*magSqr(UTilde);
        p = pStar;
    }
    else if (SNei > 0)
    {
        const scalar dS = SNei - SStar;

        this->save
        (
            facei,
            patchi,
            (SNei*rhoUNei - rhoUPhiNei + pStarStarStar*normal)
           /(rhoNei*(SNei - UvNei)),
            Uf_
        );
        phi = SStar*(SNei - UvNei)/dS;
        rhoPhi = phi*rhoNei + phip;
        rhoUPhi =
            (
                SStar*(SNei*rhoUNei - rhoUPhiNei) + SNei*pStarStarStar*normal
            )/dS
          + phip*UTilde;
        rhoEPhi =
            SStar*(SNei*rhoENei - rhoEPhiNei + SNei*pStar)/dS
          + 0.5*phip*magSqr(UTilde);
        p = pStar;
    }
    else
    {
        this->save(facei, patchi, UNei, Uf_);
        phi = UvNei;
        rhoPhi = rhoNei*UvNei;
        rhoUPhi = rhoUPhiNei;
        rhoEPhi = rhoEPhiNei;
        p = pNei;
    }

    phi *= magSf;
    rhoPhi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    rhoEPhi += meshPhi(facei, patchi)*p;
}


void Foam::fluxSchemes::HLLCP::calculateFluxes
(
    const scalarList& alphasOwn, const scalarList& alphasNei,
    const scalarList& rhosOwn, const scalarList& rhosNei,
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const scalar& cOwn, const scalar& cNei,
    const vector& Sf,
    scalar& phi,
    scalarList& alphaPhis,
    scalarList& alphaRhoPhis,
    vector& rhoUPhi,
    scalar& rhoEPhi,
    const label facei, const label patchi
)
{
    NotImplemented;

    scalar magSf = mag(Sf);
    vector normal = Sf/magSf;

    scalar EOwn = eOwn + 0.5*magSqr(UOwn);
    scalar ENei = eNei + 0.5*magSqr(UNei);

    const scalar vMesh(meshPhi(facei, patchi)/magSf);
    scalar UvOwn((UOwn & normal) - vMesh);
    scalar UvNei((UNei & normal) - vMesh);

    scalar MaOwn(UvOwn/cOwn);
    scalar MaNei(UvNei/cNei);

    scalar wOwn(sqrt(rhoOwn)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar wNei(sqrt(rhoNei)/(sqrt(rhoOwn) + sqrt(rhoNei)));
    scalar MaTilde(wOwn*MaOwn + wNei*MaNei);
    vector UTilde(wOwn*UOwn + wNei*UNei);
    scalar UvTilde(wOwn*UvOwn + wNei*UvNei);
    scalar cTilde(wOwn*cOwn + wNei*cNei);

    scalar SOwn(min(UvOwn - cOwn, UvTilde - cTilde));
    scalar SNei(max(UvNei + cNei, UvTilde + cTilde));

    scalar aOwn(rhoOwn*(SOwn - UvOwn));
    scalar aNei(rhoNei*(SNei - UvNei));

    scalar SStar((aNei*UvNei - aOwn*UvOwn + pOwn - pNei)/(aNei - aOwn));
    scalar pStar
    (
        (
            aNei*pOwn - aOwn*pNei - aOwn*aNei*(UvOwn - UvNei)
        )/(aNei - aOwn)
    );

    scalar theta(min(max(mag(MaOwn), mag(MaNei)), 1.0));
    scalar pAvg(0.5*(pOwn + pNei));
    scalar pStarStar(pStar*theta + (1.0 - theta)*pAvg);

    scalar f(f_()[facei]);
    scalar pStarStarStar(f*pStarStar + (1.0 - f)*pStar);

    scalar phip =
        (f - 1.0)
       *SOwn*SNei/(SNei - SOwn)
       /(1.0 + mag(MaTilde))
       *(pNei - pOwn)/sqr(cTilde);

    this->save(facei, patchi, SOwn, SOwn_);
    this->save(facei, patchi, SNei, SNei_);
    this->save(facei, patchi, SStar, SStar_);
    this->save(facei, patchi, pStar, pStar_);
    this->save(facei, patchi, phip, phip_);
    this->save(facei, patchi, UTilde, UTilde_);
    this->save(facei, patchi, UvOwn, UvOwn_);
    this->save(facei, patchi, UvNei, UvNei_);

    // Owner values
    const vector rhoUOwn = rhoOwn*UOwn;
    const scalar rhoEOwn = rhoOwn*EOwn;

    const vector rhoUPhiOwn = rhoUOwn*UvOwn + pOwn*normal;
    const scalar rhoEPhiOwn = (rhoEOwn + pOwn)*UvOwn;

    // Neighbour values
    const vector rhoUNei = rhoNei*UNei;
    const scalar rhoENei = rhoNei*ENei;

    const vector rhoUPhiNei = rhoUNei*UvNei + pNei*normal;
    const scalar rhoEPhiNei = (rhoENei + pNei)*UvNei;

    scalar p;
    if (SOwn > 0)
    {
        this->save(facei, patchi, UOwn, Uf_);
        phi = UvOwn;
        rhoUPhi = rhoUPhiOwn;
        rhoEPhi = rhoEPhiOwn;
        p = pOwn;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*phi;
            alphaRhoPhis[phasei] = alphaPhis[phasei]*rhosOwn[phasei];
        }
    }
    else if (SStar > 0)
    {
        const scalar dS = SOwn - SStar;

        this->save
        (
            facei,
            patchi,
            (SOwn*rhoUOwn - rhoUPhiOwn + pStarStarStar*normal)
           /(rhoOwn*(SOwn - UvOwn)),
            Uf_
        );

        scalar f = (SOwn - UvOwn)/dS;
        phi = SStar*f;

        rhoUPhi =
            (
                SStar*(SOwn*rhoUOwn - rhoUPhiOwn) + SOwn*pStarStarStar*normal
            )/dS
          + phip*UTilde;
        rhoEPhi =
            SStar*(SOwn*rhoEOwn - rhoEPhiOwn + SOwn*pStar)/dS
          + 0.5*phip*magSqr(UTilde);
        p = pStar;

        scalar phipByRho = phip/rhoOwn;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasOwn[phasei]*(phi + phipByRho);
            alphaRhoPhis[phasei] =
                alphasOwn[phasei]*rhosOwn[phasei]
               *(phi + alphasOwn[phasei]*rhosOwn[phasei]*phipByRho);
        }
    }
    else if (SNei > 0)
    {
        const scalar dS = SNei - SStar;

        this->save
        (
            facei,
            patchi,
            (SNei*rhoUNei - rhoUPhiNei + pStarStarStar*normal)
           /(rhoNei*(SNei - UvNei)),
            Uf_
        );

        scalar f = (SNei - UvNei)/dS;
        phi = SStar*f;

        rhoUPhi =
            (
                SStar*(SNei*rhoUNei - rhoUPhiNei) + SNei*pStarStarStar*normal
            )/dS
          + phip*UTilde;
        rhoEPhi =
            SStar*(SNei*rhoENei - rhoEPhiNei + SNei*pStar)/dS
          + 0.5*phip*magSqr(UTilde);
        p = pStar;

        scalar phipByRho = phip/rhoNei;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] =
                alphasNei[phasei]*(phi + phipByRho);
            alphaRhoPhis[phasei] =
                alphasNei[phasei]*rhosNei[phasei]
               *(phi + alphasNei[phasei]*rhosNei[phasei]*phipByRho);
        }
    }
    else
    {
        this->save(facei, patchi, UNei, Uf_);
        phi = UvNei;
        rhoUPhi = rhoUPhiNei;
        rhoEPhi = rhoEPhiNei;
        p = pNei;

        forAll(alphasOwn, phasei)
        {
            alphaPhis[phasei] = alphasNei[phasei]*phi;
            alphaRhoPhis[phasei] = alphaPhis[phasei]*rhosNei[phasei];
        }
    }

    phi *= magSf;
    rhoUPhi *= magSf;
    rhoEPhi *= magSf;
    rhoEPhi += meshPhi(facei, patchi)*p;
    forAll(alphasOwn, phasei)
    {
        alphaPhis[phasei] *= magSf;
        alphaRhoPhis[phasei] *= magSf;
    }
}


Foam::scalar Foam::fluxSchemes::HLLCP::energyFlux
(
    const scalar& rhoOwn, const scalar& rhoNei,
    const vector& UOwn, const vector& UNei,
    const scalar& eOwn, const scalar& eNei,
    const scalar& pOwn, const scalar& pNei,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SNei = getValue(facei, patchi, SNei_());
    scalar SStar = getValue(facei, patchi, SStar_());
    scalar pStar = getValue(facei, patchi, pStar_());
    scalar phip = getValue(facei, patchi, phip_());
    vector UTilde = getValue(facei, patchi, UTilde_());
    scalar UvOwn = getValue(facei, patchi, UvOwn_());
    scalar UvNei = getValue(facei, patchi, UvNei_());
    scalar magSf = mag(getValue(facei, patchi, mesh_.Sf()));

    // Owner values
    const scalar rhoEOwn = rhoOwn*(eOwn + 0.5*magSqr(UOwn));
    const scalar rhoEPhiOwn = (rhoEOwn + pOwn)*UvOwn;

    // Neighbour values
    const scalar rhoENei = rhoNei*(eNei + 0.5*magSqr(UNei));
    const scalar rhoEPhiNei = (rhoENei + pNei)*UvNei;
    scalar rhoEPhi, p;
    if (SOwn > 0)
    {
        rhoEPhi = rhoEPhiOwn;
        p = pOwn;
    }
    else if (SStar > 0)
    {
        const scalar dS = SOwn - SStar;
        rhoEPhi =
            SStar*(SOwn*rhoEOwn - rhoEPhiOwn + SOwn*pStar)/dS
          + 0.5*phip*magSqr(UTilde);
        p = pStar;
    }
    else if (SNei > 0)
    {
        const scalar dS = SNei - SStar;
        rhoEPhi =
            SStar*(SNei*rhoENei - rhoEPhiNei + SNei*pStar)/dS
          + 0.5*phip*magSqr(UTilde);
        p = pStar;
    }
    else
    {
        rhoEPhi = rhoEPhiNei;
        p = pNei;
    }

    return rhoEPhi*magSf + meshPhi(facei, patchi)*p;
}


Foam::scalar Foam::fluxSchemes::HLLCP::interpolate
(
    const scalar& fOwn, const scalar& fNei,
    const label facei, const label patchi
) const
{
    scalar SOwn = getValue(facei, patchi, SOwn_());
    scalar SStar = getValue(facei, patchi, SStar_());

    if (SOwn > 0 || SStar > 0)
    {
        return fOwn;
    }
    else
    {
        return fNei;
    }
}

// ************************************************************************* //
