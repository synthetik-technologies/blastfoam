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
    ),
    speedOfSound_
    (
        IOobject
        (
            phasePropertyName("speedOfSound", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimVelocity, 0.0)
    )
{}


void Foam::fluidBlastThermo::initializeFields()
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
    if (nPhases <= 1)
    {
        return basicBlastThermo::New<fluidBlastThermo>
        (
            mesh,
            dict.optionalSubDict("mixture"),
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

void Foam::fluidBlastThermo::correct()
{
    const scalarField& rhoCells = this->rho_;
    scalarField& eCells = this->e_;

    scalarField& pCells = this->p_.primitiveFieldRef();
    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& muCells = this->mu_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();
    scalarField& cCells = this->speedOfSound_.primitiveFieldRef();

    forAll(rhoCells, celli)
    {
        const blastFluidMixture& thermoMixture = this->cellThermoMixture(celli);

        // Update temperature
        TCells[celli] =
            thermoMixture.TRhoE(TCells[celli], rhoCells[celli], eCells[celli]);
        if (TCells[celli] < this->TLow_)
        {
            eCells[celli] =
                thermoMixture.HE(rhoCells[celli], eCells[celli], this->TLow_);
            TCells[celli] = TLow_;
        }

        const scalar r(rhoCells[celli]);
        const scalar e(eCells[celli]);
        const scalar T(TCells[celli]);

        pCells[celli] = thermoMixture.pRhoT(r, e, T);
        CpCells[celli] = thermoMixture.Cp(r, e, T);
        CvCells[celli] = thermoMixture.Cv(r, e, T);
        muCells[celli] = thermoMixture.mu(r, e, T);
        alphaCells[celli] = thermoMixture.kappa(r, e, T)/CpCells[celli];
        cCells[celli] = thermoMixture.speedOfSound(pCells[celli], r, e, T);
    }

    this->T_.correctBoundaryConditions();
    this->he().correctBoundaryConditions();
    this->p_.correctBoundaryConditions();


    const volScalarField::Boundary& rhoBf = this->rho_.boundaryField();
    const volScalarField::Boundary& heBf = this->he().boundaryField();
    const volScalarField::Boundary& pBf = this->p_.boundaryField();
    const volScalarField::Boundary& TBf = this->T_.boundaryField();

    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();
    volScalarField::Boundary& cBf = this->speedOfSound_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& phe = heBf[patchi];
        const fvPatchScalarField& pp = pBf[patchi];

        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];
        fvPatchScalarField& pc = cBf[patchi];

        forAll(pT, facei)
        {
            const blastFluidMixture& thermoMixture =
                this->patchFaceThermoMixture(patchi, facei);
            const scalar r(prho[facei]);
            const scalar e(phe[facei]);
            const scalar T(pT[facei]);

            pCp[facei] = thermoMixture.Cp(r, e, T);
            pCv[facei] = thermoMixture.Cv(r, e, T);
            pmu[facei] = thermoMixture.mu(r, e, T);
            palpha[facei] = thermoMixture.kappa(r, e, T)/pCp[facei];
            pc[facei] = thermoMixture.speedOfSound(pp[facei], r, e, T);
        }
    }
}



void Foam::fluidBlastThermo::updateRho()
{
    updateRho(p_);
}


Foam::volScalarField& Foam::fluidBlastThermo::p()
{
    return p_;
}


const Foam::volScalarField& Foam::fluidBlastThermo::p() const
{
    return p_;
}


const Foam::volScalarField& Foam::fluidBlastThermo::speedOfSound() const
{
    return speedOfSound_;
}


Foam::volScalarField& Foam::fluidBlastThermo::speedOfSound()
{
    return speedOfSound_;
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
