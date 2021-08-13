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

#include "singlePhaseSolidBlastThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(singlePhaseSolidBlastThermo, 0);
    addToRunTimeSelectionTable
    (
        solidBlastThermo,
        singlePhaseSolidBlastThermo,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::singlePhaseSolidBlastThermo::singlePhaseSolidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    solidBlastThermo(mesh, dict, phaseName),
    thermoPtr_
    (
        phaseSolidBlastThermo::New
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
    initializeModels();
}

void Foam::singlePhaseSolidBlastThermo::initializeModels()
{
    thermoPtr_->initializeModels();
    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        volScalarField e(thermoPtr_->calce());
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

Foam::singlePhaseSolidBlastThermo::~singlePhaseSolidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::singlePhaseSolidBlastThermo::postUpdate()
{
    thermoPtr_->postUpdate();
}


void Foam::singlePhaseSolidBlastThermo::solve()
{
    thermoPtr_->solve();
}


void Foam::singlePhaseSolidBlastThermo::update()
{

    thermoPtr_->update();
}


void Foam::singlePhaseSolidBlastThermo::updateRho()
{
    thermoPtr_->updateRho();
}


void Foam::singlePhaseSolidBlastThermo::correct()
{
    thermoPtr_->correct();

    updateRho();

    this->T_ = thermoPtr_->THE();
    this->T_.correctBoundaryConditions();
}


Foam::tmp<Foam::volScalarField>
Foam::singlePhaseSolidBlastThermo::ESource() const
{
    return thermoPtr_->ESource();
}


Foam::tmp<Foam::volVectorField>
Foam::singlePhaseSolidBlastThermo::Kappa() const
{
    return thermoPtr_->Kappa();
}


Foam::tmp<Foam::vectorField>
Foam::singlePhaseSolidBlastThermo::Kappa(const label patchi) const
{
    return thermoPtr_->Kappa(patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::he
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return thermoPtr_->he(p, T);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::he
(
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->he(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::he
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->he(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::hs() const
{
    return thermoPtr_->hs();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::hs
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return thermoPtr_->hs(p, T);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::hs
(
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->hs(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::hs
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->hs(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::ha() const
{
    return thermoPtr_->ha();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::ha
(
    const volScalarField& p,
    const volScalarField& T
) const
{
    return thermoPtr_->ha(p, T);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::ha
(
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->ha(T, cells);
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::ha
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->ha(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::hc() const
{
    return thermoPtr_->hc();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::flameT() const
{
    return thermoPtr_->flameT();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::THE() const
{
    return thermoPtr_->THE();
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::THE
(
    const volScalarField& he,
    const volScalarField& p,
    const volScalarField& T0
) const
{
    return thermoPtr_->THE(he, p, T0);
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseSolidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const labelList& cells
) const
{
    return thermoPtr_->THE(he, T, cells);
}


Foam::tmp<Foam::scalarField>
Foam::singlePhaseSolidBlastThermo::THE
(
    const scalarField& he,
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->THE(he, T, patchi);
}


Foam::scalar Foam::singlePhaseSolidBlastThermo::THEi
(
    const scalar he,
    const scalar T,
    const label celli
) const
{
    return thermoPtr_->THEi(he, T, celli);
}

Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::Cp() const
{
    return thermoPtr_->Cp();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::Cp
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->Cp(T, patchi);
}


Foam::scalar Foam::singlePhaseSolidBlastThermo::Cpi(const label celli) const
{
    return thermoPtr_->Cpi(celli);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::Cv() const
{
    return thermoPtr_->Cv();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::Cv
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->Cv(T, patchi);
}


Foam::scalar Foam::singlePhaseSolidBlastThermo::Cvi(const label celli) const
{
    return thermoPtr_->Cvi(celli);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::Cpv() const
{
    return thermoPtr_->Cpv();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::Cpv
(
    const scalarField& T,
    const label patchi
) const
{
    return thermoPtr_->Cpv(T, patchi);
}


Foam::tmp<Foam::volScalarField> Foam::singlePhaseSolidBlastThermo::W() const
{
    return thermoPtr_->W();
}


Foam::tmp<Foam::scalarField> Foam::singlePhaseSolidBlastThermo::W
(
    const label patchi
) const
{
    return thermoPtr_->W(patchi);
}


Foam::scalar Foam::singlePhaseSolidBlastThermo::Wi(const label celli) const
{
    return thermoPtr_->Wi(celli);
}

// ************************************************************************* //
