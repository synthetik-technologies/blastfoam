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

#include "twoPhaseFluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(twoPhaseFluidBlastThermo, 0);
    addToRunTimeSelectionTable
    (
        fluidBlastThermo,
        twoPhaseFluidBlastThermo,
        phase
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseFluidBlastThermo::twoPhaseFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidBlastThermo(mesh, dict, phaseName),
    phases_(dict.lookup("phases")),
    volumeFraction_
    (
        IOobject
        (
            IOobject::groupName("alpha", phases_[0]),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    rho1_
    (
        IOobject
        (
            IOobject::groupName("rho", phases_[0]),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    thermo1_
    (
        phaseFluidBlastThermo::New
        (
            mesh,
            dict.subDict(phases_[0]),
            phases_[0],
            phaseName
        )
    ),
    rho2_
    (
        IOobject
        (
            IOobject::groupName("rho", phases_[1]),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    thermo2_
    (
        phaseFluidBlastThermo::New
        (
            mesh,
            dict.subDict(phases_[1]),
            phases_[1],
            phaseName
        )
    )
{
    //- Force reading of residual values
    thermo1_->read(dict.subDict(phases_[0]));
    thermo2_->read(dict.subDict(phases_[1]));

    this->residualAlpha_ =
        max(thermo1_->residualAlpha(), thermo2_->residualAlpha());
    this->residualRho_ = max(thermo1_->residualRho(), thermo2_->residualRho());

    // Update total density
    this->rho_ = volumeFraction_*rho1_ + (1.0 - volumeFraction_)*rho2_;

    mu_ =
        volumeFraction_*thermo1_->calcMu()
      + (1.0 - volumeFraction_)*thermo2_->calcMu();
}

void Foam::twoPhaseFluidBlastThermo::initializeModels()
{
    thermo1_->initializeModels();
    thermo2_->initializeModels();

    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        volScalarField e
        (
            thermo1_->calce(this->p())*volumeFraction_
          + thermo2_->calce(this->p())*(1.0 - volumeFraction_)
        );
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

Foam::twoPhaseFluidBlastThermo::~twoPhaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::twoPhaseFluidBlastThermo::read(const dictionary& dict)
{
    thermo1_->read(dict.subDict(thermo1_->name()));
    thermo2_->read(dict.subDict(thermo2_->name()));
}


void Foam::twoPhaseFluidBlastThermo::postUpdate()
{
    thermo1_->postUpdate();
    thermo2_->postUpdate();
}


void Foam::twoPhaseFluidBlastThermo::solve()
{
    thermo1_->solve();
    thermo2_->solve();
}


void Foam::twoPhaseFluidBlastThermo::update()
{

    thermo1_->update();
    thermo2_->update();
}


void Foam::twoPhaseFluidBlastThermo::updateRho(const volScalarField& p)
{
    thermo1_->updateRho(p);
    thermo2_->updateRho(p);
    this->rho_ =
        volumeFraction_*thermo1_->rho()
      + (1.0 - volumeFraction_)*thermo2_->rho();
}


void Foam::twoPhaseFluidBlastThermo::correct()
{
    this->T_ =
        volumeFraction_*thermo1_->THE()
      + (1.0 - volumeFraction_)*thermo2_->THE();
    this->T_.correctBoundaryConditions();

    {
        volScalarField alphaXi1
        (
            volumeFraction_/(thermo1_->Gamma() - 1.0)
        );
        volScalarField alphaXi2
        (
            (1.0 - volumeFraction_)/(thermo2_->Gamma() - 1.0)
        );

        this->p() =
            (
                alphaXi1*thermo1_->pRhoT()*pos(alphaXi1 - small)
              + alphaXi2*thermo2_->pRhoT()*pos(alphaXi2 - small)
            )/(alphaXi1 + alphaXi2);
        this->p().correctBoundaryConditions();
    }

    // Update transport coefficients
    mu_ =
        volumeFraction_*thermo1_->calcMu()
      + (1.0 - volumeFraction_)*thermo2_->calcMu();
    this->alpha_ = this->kappa()/this->Cp();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::ESource() const
{
    return
        volumeFraction_*thermo1_->ESource()
      + (1.0 - volumeFraction_)*thermo2_->ESource();
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::calce(const volScalarField& p) const
{
    return
        volumeFraction_*thermo1_->calce(p)
      + (1.0 - volumeFraction_)*thermo2_->calce(p);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::pRhoT() const
{

    volScalarField alphaXi1
    (
        volumeFraction_/(thermo1_->Gamma() - 1.0)
    );
    volScalarField alphaXi2
    (
        (1.0 - volumeFraction_)/(thermo2_->Gamma() - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->pRhoT()*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->pRhoT()*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellpRhoT(const label celli) const
{
    scalar alphaXi1
    (
        volumeFraction_[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        (1.0 - volumeFraction_[celli])/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->cellpRhoT(celli)*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->cellpRhoT(celli)*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpdRho(const label celli) const
{
    scalar alphaXi1
    (
        volumeFraction_[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        (1.0 - volumeFraction_[celli])/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpdRho(celli)*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->celldpdRho(celli)*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::celldpde(const label celli) const
{
    scalar alphaXi1
    (
        volumeFraction_[celli]/(thermo1_->cellGamma(celli) - 1.0)
    );
    scalar alphaXi2
    (
        (1.0 - volumeFraction_[celli])/(thermo2_->cellGamma(celli) - 1.0)
    );

    return
        (
            alphaXi1*thermo1_->celldpde(celli)*pos(alphaXi1 - small)
          + alphaXi2*thermo2_->celldpde(celli)*pos(alphaXi2 - small)
        )/(alphaXi1 + alphaXi2);
}


Foam::tmp<Foam::volScalarField>
Foam::twoPhaseFluidBlastThermo::speedOfSound() const
{
    volScalarField alphaXi1(volumeFraction_/(thermo1_->Gamma() - 1.0));
    volScalarField alphaXi2((1.0 - volumeFraction_)/(thermo2_->Gamma() - 1.0));

    volScalarField cSqr
    (
        (
            alphaXi1*rho1_*sqr(thermo1_->speedOfSound())
          + alphaXi2*rho2_*sqr(thermo2_->speedOfSound())
        )
        /(alphaXi1 + alphaXi2)/this->rho_
    );
    cSqr.max(small);
    return sqrt(cSqr);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::speedOfSound(const label patchi) const
{
    scalarField alphaXi1
    (
        volumeFraction_.boundaryField()[patchi]/(thermo1_->Gamma(patchi) - 1.0)
    );
    scalarField alphaXi2
    (
        (1.0 - volumeFraction_.boundaryField()[patchi])
       /(thermo2_->Gamma(patchi) - 1.0)
    );

    scalarField cSqr
    (
        (
            alphaXi1*rho1_.boundaryField()[patchi]
           *sqr(thermo1_->speedOfSound(patchi))
          + alphaXi2*rho2_.boundaryField()[patchi]
           *sqr(thermo2_->speedOfSound(patchi))
        )
        /(this->rho_.boundaryField()[patchi]*(alphaXi1 + alphaXi2)));
    return sqrt(max(cSqr, small));
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::Gamma() const
{
    return
        volumeFraction_*thermo1_->Gamma()
      + (1.0 - volumeFraction_)*thermo2_->Gamma();
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::Gamma(const label patchi) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->Gamma(patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->Gamma(patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellGamma(const label celli) const
{
    return
        volumeFraction_[celli]*thermo1_->cellGamma(celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->cellGamma(celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return
        volumeFraction_*thermo1_->he(p, T)
      + (1.0 - volumeFraction_)*thermo2_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    scalarField vf(UIndirectList<scalar>(volumeFraction_(), cells)());
    return
        vf*thermo1_->he(T, cells)
      + (1.0 - vf)*thermo2_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->he(T, patchi)
      + (1.0 - vf)*thermo2_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hs() const
{
    return
        volumeFraction_*thermo1_->hs()
      + (1.0 - volumeFraction_)*thermo2_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return
        volumeFraction_*thermo1_->hs(p, T)
      + (1.0 - volumeFraction_)*thermo2_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    scalarField vf(UIndirectList<scalar>(volumeFraction_(), cells)());
    return
        vf*thermo1_->hs(T, cells)
      + (1.0 - vf)*thermo2_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->hs(T, patchi)
      + (1.0 - vf)*thermo2_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::ha() const
{
    return
        volumeFraction_*thermo1_->ha()
      + (1.0 - volumeFraction_)*thermo2_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return
        volumeFraction_*thermo1_->ha(p, T)
      + (1.0 - volumeFraction_)*thermo2_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    scalarField vf(UIndirectList<scalar>(volumeFraction_(), cells)());
    return
        vf*thermo1_->ha(T, cells)
      + (1.0 - vf)*thermo2_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->ha(T, patchi)
      + (1.0 - vf)*thermo2_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::hc() const
{
    return
        volumeFraction_*thermo1_->hc()
      + (1.0 - volumeFraction_)*thermo2_->hc();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::flameT() const
{
    return
        volumeFraction_*thermo1_->flameT()
      + (1.0 - volumeFraction_)*thermo2_->flameT();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::THE() const
{
    return
        volumeFraction_*thermo1_->THE()
      + (1.0 - volumeFraction_)*thermo2_->THE();
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return
        volumeFraction_*thermo1_->THE(he, p, T0)
      + (1.0 - volumeFraction_)*thermo2_->THE(he, p, T0);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    scalarField vf(UIndirectList<scalar>(volumeFraction_(), cells)());
    return
        vf*thermo1_->THE(he, T, cells)
      + (1.0 - vf)*thermo2_->THE(he, T, cells);
}


Foam::tmp<Foam::scalarField>
Foam::twoPhaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return
        volumeFraction_.boundaryField()[patchi]
       *thermo1_->THE(he, T, patchi)
      + (1.0 - volumeFraction_.boundaryField()[patchi])
       *thermo2_->THE(he, T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellTHE
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return
        volumeFraction_[celli]*thermo1_->cellTHE(he, T, celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->cellTHE(he, T, celli);
}

Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::Cp() const
{
    return
        volumeFraction_*thermo1_->Cp()
      + (1.0 - volumeFraction_)*thermo2_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->Cp(T, patchi)
      + (1.0 - vf)*thermo2_->Cp(T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCp
(
    const scalar T,
    const label celli
) const
{
    return
        volumeFraction_[celli]*thermo1_->cellCp(T, celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->cellCp(T, celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::Cv() const
{
    return
        volumeFraction_*thermo1_->Cv()
      + (1.0 - volumeFraction_)*thermo2_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->Cv(T, patchi)
      + (1.0 - vf)*thermo2_->Cv(T, patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellCv
(
    const scalar T,
    const label celli
) const
{
    return
        volumeFraction_[celli]*thermo1_->cellCv(T, celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->cellCv(T, celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::Cpv() const
{
    return
        volumeFraction_*thermo1_->Cpv()
      + (1.0 - volumeFraction_)*thermo2_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->Cpv(T, patchi)
      + (1.0 - vf)*thermo2_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::W() const
{
    return
        volumeFraction_*thermo1_->W()
      + (1.0 - volumeFraction_)*thermo2_->W();
}


Foam::tmp<Foam::scalarField> Foam::twoPhaseFluidBlastThermo::W
(
    const label patchi
) const
{
    const scalarField& vf(volumeFraction_.boundaryField()[patchi]);
    return
        vf*thermo1_->W(patchi)
      + (1.0 - vf)*thermo2_->W(patchi);
}


Foam::scalar Foam::twoPhaseFluidBlastThermo::cellW(const label celli) const
{
    return
        volumeFraction_[celli]*thermo1_->cellW(celli)
      + (1.0 - volumeFraction_[celli])*thermo2_->cellW(celli);
}


Foam::tmp<Foam::volScalarField> Foam::twoPhaseFluidBlastThermo::calcMu() const
{
    return
        volumeFraction_*thermo1_->calcMu()
      + (1.0 - volumeFraction_)*thermo2_->calcMu();
}

// ************************************************************************* //
