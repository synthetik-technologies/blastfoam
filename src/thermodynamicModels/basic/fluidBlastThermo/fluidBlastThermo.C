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

#include "fluidBlastThermo.H"
#include "basicBlastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidBlastThermo, 0);
    defineRunTimeSelectionTable(fluidBlastThermo, dictionary);
    defineRunTimeSelectionTable(fluidBlastThermo, phase);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidBlastThermo::fluidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word&
)
:
    blastThermo(mesh, dict, phaseName),
    p_
    (
        basicBlastThermo::lookupOrConstruct
        (
            mesh,
            phasePropertyName("p", phaseName),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE,
            dimPressure,
            true // pressure is always allowed to read from mixture field
        )
    ),
    mu_
    (
        IOobject
        (
            phasePropertyName("thermo:mu", phaseName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimensionSet(1, -1, -1, 0, 0), 0.0)
    )
{}


void Foam::fluidBlastThermo::initializeModels()
{
    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        volScalarField e(this->calce(p_));
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


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fluidBlastThermo> Foam::fluidBlastThermo::New
(
    const label nPhases,
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    if (nPhases == 1)
    {
        return basicBlastThermo::New<fluidBlastThermo>
        (
            mesh,
            dict,
            phaseName,
            phaseName
        );
    }

    word fluidType("twoPhaseFluid");
    if (nPhases > 2)
    {
        fluidType = "multiphaseFluid";
    }

    phaseConstructorTable::iterator cstrIter =
        phaseConstructorTablePtr_->find(fluidType);

    if (cstrIter == phaseConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown number of fluids " << endl
            << phaseConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return cstrIter()(mesh, dict, phaseName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidBlastThermo::~fluidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidBlastThermo::updateRho()
{
    updateRho(p_);
}


void Foam::fluidBlastThermo::correct()
{
    this->T_ = this->THE();
    this->T_.correctBoundaryConditions();

    this->p() = this->pRhoT();
    this->p().correctBoundaryConditions();

    this->mu_ = this->calcMu();
    this->alpha_ = this->kappa()/this->Cp();
}


Foam::volScalarField& Foam::fluidBlastThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::fluidBlastThermo::p() const
{
    return p_;
}


Foam::tmp<Foam::volScalarField> Foam::fluidBlastThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField> Foam::fluidBlastThermo::mu
(
    const label patchi
) const
{
    return mu_.boundaryField()[patchi];
}

Foam::scalar Foam::fluidBlastThermo::cellnu
(
    const label celli
) const
{
    return mu_[celli]/cellrho(celli);
}
// ************************************************************************* //
