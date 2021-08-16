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

#include "multiphaseFluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiphaseFluidBlastThermo, 0);
    addToRunTimeSelectionTable
    (
        fluidBlastThermo,
        multiphaseFluidBlastThermo,
        phase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiphaseFluidBlastThermo::multiphaseFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidBlastThermo(mesh, dict, phaseName),
    phases_(dict.lookup("phases")),
    volumeFractions_(phases_.size()),
    rhos_(phases_.size()),
    thermos_(phases_.size()),
    alphaRhos_(phases_.size()),
    residualAlpha_(dimless, 0.0),
    sumVfPtr_(nullptr)
{
    volScalarField sumAlpha
    (
        IOobject
        (
            "sumAlpha",
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        0.0
    );
    forAll(phases_, phasei)
    {
        word phaseIName = phases_[phasei];
        volumeFractions_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("alpha", phases_[phasei]),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        rhos_.set
        (
            phasei,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("rho", phases_[phasei]),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
        thermos_.set
        (
            phasei,
            phaseFluidBlastThermo::New
            (
                mesh,
                dict.subDict(phaseIName),
                phaseIName,
                phaseName
            ).ptr()
        );
        thermos_[phasei].read(dict.subDict(phaseIName));

        residualAlpha_ = max(residualAlpha_, thermos_[phasei].residualAlpha());

        mu_ += volumeFractions_[phasei]*thermos_[phasei].mu();

        rho_ += volumeFractions_[phasei]*rhos_[phasei];
        sumAlpha += volumeFractions_[phasei];
    }
    rho_ /= max(sumAlpha, residualAlpha_);
}


void Foam::multiphaseFluidBlastThermo::initializeModels()
{
    forAll(thermos_, phasei)
    {
        thermos_[phasei].initializeModels();
    }

    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        volScalarField e(volumeFractions_[0]*thermos_[0].calce(this->p()));
        for (label phasei = 1; phasei < thermos_.size(); phasei++)
        {
            e += volumeFractions_[phasei]*thermos_[phasei].calce(this->p());
        }
        normalise(e);

        e_ = e;

        //- Force fixed boundaries to be updates
        forAll(e_.boundaryField(), patchi)
        {
            forAll(e_.boundaryField()[patchi], facei)
            {
                e_.boundaryFieldRef()[patchi][facei] =
                    e.boundaryField()[patchi][facei];
            }
        }
    }

    correct();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiphaseFluidBlastThermo::~multiphaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiphaseFluidBlastThermo::read(const dictionary& dict)
{
    forAll(thermos_, phasei)
    {
        thermos_[phasei].read(dict.subDict(thermos_[phasei].name()));
    }
}


void Foam::multiphaseFluidBlastThermo::setTotalVolumeFractionPtr
(
    const volScalarField& vf
)
{
    sumVfPtr_.set(&vf);
}


void Foam::multiphaseFluidBlastThermo::update()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].update();
    }
}


void Foam::multiphaseFluidBlastThermo::solve()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].solve();
    }
}


void Foam::multiphaseFluidBlastThermo::postUpdate()
{
    forAll(phases_, phasei)
    {
        thermos_[phasei].postUpdate();
    }
}


void Foam::multiphaseFluidBlastThermo::updateRho(const volScalarField& p)
{
    thermos_[0].updateRho(p);
    this->rho_ = volumeFractions_[0]*thermos_[0].rho();
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        thermos_[phasei].updateRho(p);
        rho_ += volumeFractions_[phasei]*thermos_[phasei].rho();
    }
    normalise(rho_);
}


void Foam::multiphaseFluidBlastThermo::correct()
{
    volScalarField rGamma(volumeFractions_[0]/(thermos_[0].Gamma() - 1.0));
    volScalarField pByGamma
    (
        rGamma*thermos_[0].pRhoT()*pos(volumeFractions_[0] - residualAlpha_)
    );

    T_ = volumeFractions_[0]*thermos_[0].THE();
    mu_ = volumeFractions_[0]*thermos_[0].mu();
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        const volScalarField& vf(volumeFractions_[phasei]);
        volScalarField alphaByGamma(vf/(thermos_[phasei].Gamma() - 1.0));
        rGamma += alphaByGamma;
        pByGamma +=
            alphaByGamma*thermos_[phasei].pRhoT()*pos(vf - residualAlpha_);

        T_ += vf*thermos_[phasei].THE();
        mu_ += vf*thermos_[phasei].mu();
    }
    p_ = pByGamma/max(rGamma, residualAlpha_);
    p_.correctBoundaryConditions();

    normalise(T_);
    T_.correctBoundaryConditions();

    normalise(mu_);
    this->alpha_ = this->kappa()/this->Cp();
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidBlastThermo::ESource() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].ESource());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ESource();
    }
    return tmpF;
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::calce(const volScalarField& p) const
{
    tmp<volScalarField> e(volumeFractions_[0]*thermos_[0].calce(p));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        e.ref() += volumeFractions_[phasei]*thermos_[phasei].calce(p);
    }
    return normalise(e);
}


Foam::tmp<Foam::volScalarField> Foam::multiphaseFluidBlastThermo::pRhoT() const
{
    volScalarField rGamma(volumeFractions_[0]/(thermos_[0].Gamma() - 1.0));
    volScalarField pByGamma
    (
        rGamma*thermos_[0].pRhoT()*pos(volumeFractions_[0] - residualAlpha_)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        const volScalarField& vf(volumeFractions_[phasei]);
        volScalarField alphaByGamma(vf/(thermos_[phasei].Gamma() - 1.0));
        rGamma += alphaByGamma;
        pByGamma +=
            alphaByGamma*thermos_[phasei].pRhoT()*pos(vf - residualAlpha_);
    }
    return pByGamma/max(rGamma, residualAlpha_);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::pRhoTi(const label celli) const
{
    scalar rGamma = 0.0;
    scalar pByGamma = 0.0;
    forAll(thermos_, phasei)
    {
        scalar vf(volumeFractions_[phasei][celli]);
        scalar alphaByGamma(vf/(thermos_[phasei].Gammai(celli) - 1.0));
        rGamma += alphaByGamma;
        pByGamma +=
            alphaByGamma
           *thermos_[phasei].pRhoTi(celli)*pos(vf - residualAlpha_.value());
    }
    return pByGamma/max(rGamma, residualAlpha_.value());
}


Foam::scalar Foam::multiphaseFluidBlastThermo::dpdRhoi(const label celli) const
{
    scalar rGamma = 0.0;
    scalar pByGamma = 0.0;
    forAll(thermos_, phasei)
    {
        scalar vf(volumeFractions_[phasei][celli]);
        scalar alphaByGamma(vf/(thermos_[phasei].Gammai(celli) - 1.0));
        rGamma += alphaByGamma;
        pByGamma +=
            alphaByGamma
           *thermos_[phasei].dpdRhoi(celli)*pos(vf - residualAlpha_.value());
    }
    return pByGamma/max(rGamma, residualAlpha_.value());
}


Foam::scalar Foam::multiphaseFluidBlastThermo::dpdei(const label celli) const
{
    scalar rGamma = 0.0;
    scalar pByGamma = 0.0;
    forAll(thermos_, phasei)
    {
        scalar vf(volumeFractions_[phasei][celli]);
        scalar alphaByGamma(vf/(thermos_[phasei].Gammai(celli) - 1.0));
        rGamma += alphaByGamma;
        pByGamma +=
            alphaByGamma
           *thermos_[phasei].dpdei(celli)*pos(vf - residualAlpha_.value());
    }
    return pByGamma/max(rGamma, residualAlpha_.value());
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::speedOfSound() const
{
    volScalarField alphaRhoXii
    (
        volumeFractions_[0]*rhos_[0]/(thermos_[0].Gamma() - 1.0)
    );
    volScalarField alphaXiRhoCSqr
    (
        "alphaXiRhoCSqr",
         alphaRhoXii*sqr(thermos_[0].speedOfSound())
    );
    volScalarField alphaRhoXi("alphaRhoXi", alphaRhoXii);

    forAll(thermos_, phasei)
    {
        alphaRhoXii =
            volumeFractions_[phasei]
           *rhos_[phasei]/(thermos_[phasei].Gamma() - 1.0);
        alphaRhoXi += alphaRhoXii;
        alphaXiRhoCSqr += alphaRhoXii*sqr(thermos_[phasei].speedOfSound());
    }
    alphaRhoXi.max(1e-10);
    tmp<volScalarField> cSqr(alphaXiRhoCSqr/alphaRhoXi);
    cSqr.ref().max(small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::speedOfSound(const label patchi) const
{
    scalarField alphaXiRhoCSqr(p_.boundaryField()[patchi].size(), 0.0);
    scalarField alphaRhoXi(p_.boundaryField()[patchi].size(), 0.0);

    forAll(thermos_, phasei)
    {
        scalarField alphaRhoXii
        (
            volumeFractions_[phasei].boundaryField()[patchi]
           *rhos_[phasei].boundaryField()[patchi]
           /(thermos_[phasei].Gamma(patchi) - 1.0)
        );
        alphaRhoXi += alphaRhoXii;
        alphaXiRhoCSqr +=
            alphaRhoXii
           *sqr(thermos_[phasei].speedOfSound(patchi));
    }

    scalarField cSqr(alphaXiRhoCSqr/max(alphaRhoXi, 1e-10));
    cSqr = max(cSqr, small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::Gamma() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]/thermos_[0].Gamma());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]/thermos_[phasei].Gamma();
    }
    return 1.0/normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Gamma(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       /thermos_[0].Gamma(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           /thermos_[phasei].Gamma(patchi);
    }
    return 1.0/normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::Gammai(const label celli) const
{
    scalar f(volumeFractions_[0][celli]/thermos_[0].Gammai(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]/thermos_[phasei].Gammai(celli);
    }
    return 1.0/normalise(f, celli);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].he(p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].he(p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].he(T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].he(T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].he(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].he(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::hs() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].hs());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].hs();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].hs(p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].hs(p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].hs(T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].hs(T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].hs(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].hs(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::ha() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].ha());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ha();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].ha(p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].ha(p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].ha(T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].ha(T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].ha(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].ha(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::hc() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].hc());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].hc();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::flameT() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].flameT());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].flameT();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::THE() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].THE());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].THE();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T
) const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].THE(he, p, T));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].THE(he, p, T);
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    tmp<scalarField> tmpF
    (
        UIndirectList<scalar>(volumeFractions_[0], cells)()
       *thermos_[0].THE(he, T, cells)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            UIndirectList<scalar>(volumeFractions_[phasei], cells)()
           *thermos_[phasei].THE(he, T, cells);
    }
    return normalise(tmpF, cells);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].THE(he, T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].THE(he, T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::THEi
(
    const scalar he,
    const scalar T0,
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].THEi(he, T0, celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].THEi(he, T0, celli);
    }
    return normalise(f, celli);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::Cp() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].Cp());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].Cp();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cp(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cp(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::Cpi
(
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].Cpi(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].Cpi(celli);
    }
    return normalise(f, celli);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::Cv() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].Cv());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].Cv();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cp(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cv(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::Cvi
(
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].Cvi(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].Cvi(celli);
    }
    return normalise(f, celli);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::Cpv() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].Cpv());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].Cpv();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].Cpv(T, patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].Cpv(T, patchi);
    }
    return normalise(tmpF, patchi);
}


Foam::tmp<Foam::volScalarField>
Foam::multiphaseFluidBlastThermo::W() const
{
    tmp<volScalarField> tmpF(volumeFractions_[0]*thermos_[0].W());
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() += volumeFractions_[phasei]*thermos_[phasei].W();
    }
    return normalise(tmpF);
}


Foam::tmp<Foam::scalarField>
Foam::multiphaseFluidBlastThermo::W(const label patchi) const
{
    tmp<scalarField> tmpF
    (
        volumeFractions_[0].boundaryField()[patchi]
       *thermos_[0].W(patchi)
    );
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        tmpF.ref() +=
            volumeFractions_[phasei].boundaryField()[patchi]
           *thermos_[phasei].W(patchi);
    }
    return 1.0/normalise(tmpF, patchi);
}


Foam::scalar Foam::multiphaseFluidBlastThermo::Wi
(
    const label celli
) const
{
    scalar f(volumeFractions_[0][celli]*thermos_[0].Wi(celli));
    for (label phasei = 1; phasei < thermos_.size(); phasei++)
    {
        f += volumeFractions_[phasei][celli]*thermos_[phasei].Wi(celli);
    }
    return normalise(f, celli);
}


// ************************************************************************* //
