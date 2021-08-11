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

#include "singlePhaseFluidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseFluidBlastThermo, 0);
    addToRunTimeSelectionTable
    (
        fluidBlastThermo,
        singlePhaseFluidBlastThermo,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseFluidBlastThermo::singlePhaseFluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    fluidBlastThermo(mesh, dict, phaseName),
    thermoPtr_
    (
        phaseFluidBlastThermo::New
        (
            phaseName,
            blastThermo::rho(),
            blastThermo::he(),
            T(),
            dict.optionalSubDict("mixture"),
            phaseName
        )
    )
{
    if (!rho().typeHeaderOk<volScalarField>(true))
    {
        FatalErrorInFunction
            << rho().name() << " must be read." << nl
            << abort(FatalError);
    }

    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        volScalarField e(thermoPtr_->calce(this->p()));
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
        e_.correctBoundaryConditions();
    }
    correct();
}

void Foam::singlePhaseFluidBlastThermo::initializeModels()
{
    thermoPtr_->initializeModels();
    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        volScalarField e(thermoPtr_->calce(this->p()));
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
        e_.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::singlePhaseFluidBlastThermo::~singlePhaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singlePhaseFluidBlastThermo::postUpdate()
{
    thermoPtr_->postUpdate();
}


void Foam::singlePhaseFluidBlastThermo::solve()
{
    thermoPtr_->solve();
}


void Foam::singlePhaseFluidBlastThermo::update()
{

    thermoPtr_->update();
}


void Foam::singlePhaseFluidBlastThermo::updateRho(const volScalarField& p)
{
    thermoPtr_->updateRho(p);
}


void Foam::singlePhaseFluidBlastThermo::correct()
{

    this->T_ = thermoPtr_->THE();
    this->T_.correctBoundaryConditions();

    this->p() = thermoPtr_->pRhoT();
    this->p().correctBoundaryConditions();

    thermoPtr_->correct();
}


Foam::tmp<Foam::volScalarField>
Foam::singlePhaseFluidBlastThermo::ESource() const
{
    return thermoPtr_->ESource();
}


Foam::tmp<Foam::volScalarField>
Foam::singlePhaseFluidBlastThermo::speedOfSound() const
{
    return thermoPtr_->speedOfSound();
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseFluidBlastThermo::speedOfSound(const label patchi) const
{
    return thermoPtr_->speedOfSound(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::Gamma() const
{
    return thermoPtr_->Gamma();
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseFluidBlastThermo::Gamma(const label patchi) const
{
    return thermoPtr_->Gamma(patchi);
}


Foam::scalar Foam::singlePhaseFluidBlastThermo::Gammai(const label celli) const
{
    return thermoPtr_->Gammai(celli);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return thermoPtr_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::hs() const
{
    return thermoPtr_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return thermoPtr_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::ha() const
{
    return thermoPtr_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return thermoPtr_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::hc() const
{
    return thermoPtr_->hc();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::flameT() const
{
    return thermoPtr_->flameT();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::THE() const
{
    return thermoPtr_->THE();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return thermoPtr_->THE(he, p, T0);
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->THE(he, T, cells);
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseFluidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->THE(he, T, patchi);
}


Foam::scalar Foam::singlePhaseFluidBlastThermo::THEi
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return thermoPtr_->THEi(he, T, celli);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::Cp() const
{
    return thermoPtr_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->Cp(T, patchi);
}


Foam::scalar Foam::singlePhaseFluidBlastThermo::Cpi(const label celli) const
{
    return thermoPtr_->Cpi(celli);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::Cv() const
{
    return thermoPtr_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->Cv(T, patchi);
}


Foam::scalar Foam::singlePhaseFluidBlastThermo::Cvi(const label celli) const
{
    return thermoPtr_->Cvi(celli);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::Cpv() const
{
    return thermoPtr_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseFluidBlastThermo::W() const
{
    return thermoPtr_->W();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseFluidBlastThermo::W
(
    const label patchi
) const
{
    return thermoPtr_->W(patchi);
}


Foam::scalar Foam::singlePhaseFluidBlastThermo::Wi(const label celli) const
{
    return thermoPtr_->Wi(celli);
}

// ************************************************************************* //
