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

#include "solidBlastThermo.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "coordinateSystem.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidBlastThermo, 0);
    defineRunTimeSelectionTable(solidBlastThermo, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBlastThermo::solidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    blastThermo(mesh, dict, phaseName)
{}


void Foam::solidBlastThermo::initializeFields()
{
    updateRho();
    if (!e_.typeHeaderOk<volScalarField>(true))
    {
        e_ == this->calce();
    }

    correct();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBlastThermo> Foam::solidBlastThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    return blastThermo::New<solidBlastThermo>
    (
        mesh,
        dict,
        phaseName,
        phaseName
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBlastThermo::~solidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
/*
void Foam::solidBlastThermo::correct()
{
    const scalarField& rhoCells = this->rho_;
    scalarField& eCells = this->e_;

    scalarField& TCells = this->T_.primitiveFieldRef();
    scalarField& CpCells = this->Cp_.primitiveFieldRef();
    scalarField& CvCells = this->Cv_.primitiveFieldRef();
    scalarField& alphaCells = this->alpha_.primitiveFieldRef();

    forAll(rhoCells, celli)
    {
        const blastMixture& thermoMixture = this->cellThermoMixture(celli);

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

        CpCells[celli] = thermoMixture.Cp(r, e, T);
        CvCells[celli] = thermoMixture.Cv(r, e, T);
        alphaCells[celli] = thermoMixture.kappa(r, e, T)/CvCells[celli];
    }

    this->T_.correctBoundaryConditions();
    this->he().correctBoundaryConditions();

    const volScalarField::Boundary& rhoBf = this->rho_.boundaryField();
    const volScalarField::Boundary& heBf = this->he().boundaryFieldRef();

    volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
    volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
    volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = this->alpha_.boundaryFieldRef();

    forAll(this->T_.boundaryField(), patchi)
    {
        const fvPatchScalarField& prho = rhoBf[patchi];
        const fvPatchScalarField& pT = TBf[patchi];
        const fvPatchScalarField& phe = heBf[patchi];

        fvPatchScalarField& pCp = CpBf[patchi];
        fvPatchScalarField& pCv = CvBf[patchi];
        fvPatchScalarField& palpha = alphaBf[patchi];

        forAll(pT, facei)
        {
            const blastMixture& thermoMixture =
                this->patchFaceThermoMixture(patchi, facei);
            const scalar r(prho[facei]);
            const scalar e(phe[facei]);
            const scalar T(pT[facei]);

            pCp[facei] = thermoMixture.Cp(r, e, T);
            pCv[facei] = thermoMixture.Cv(r, e, T);
            palpha[facei] = thermoMixture.kappa(r, e, T)/pCp[facei];
        }
    }
}*/


Foam::tmp<Foam::volScalarField> Foam::solidBlastThermo::nu() const
{
    return volScalarField::New
    (
        "nu",
        this->T_.mesh(),
        dimViscosity
    );
}


Foam::tmp<Foam::volSymmTensorField> Foam::solidBlastThermo::KappaLocal() const
{
    const fvMesh& mesh = this->T_.mesh();

    const coordinateSystem coordinates
    (
        coordinateSystem::New(mesh, this->properties())
    );

    const tmp<volVectorField> tKappa(Kappa());
    const volVectorField& Kappa = tKappa();

    tmp<volSymmTensorField> tKappaLocal
    (
        volSymmTensorField::New
        (
            "KappaLocal",
            mesh,
            dimensionedSymmTensor(Kappa.dimensions(), Zero)
        )
    );
    volSymmTensorField& KappaLocal = tKappaLocal.ref();

    KappaLocal.primitiveFieldRef() =
        coordinates.R(mesh.C()).transformVector(Kappa);

    forAll(KappaLocal.boundaryField(), patchi)
    {
        KappaLocal.boundaryFieldRef()[patchi] =
            coordinates.R(mesh.boundary()[patchi].Cf())
           .transformVector(Kappa.boundaryField()[patchi]);
    }

    return tKappaLocal;
}


Foam::tmp<Foam::symmTensorField> Foam::solidBlastThermo::KappaLocal
(
    const label patchi
) const
{
    const fvMesh& mesh = this->T_.mesh();

    const coordinateSystem coordinates
    (
        coordinateSystem::New(mesh, this->properties())
    );

    return
        coordinates.R(mesh.boundary()[patchi].Cf())
       .transformVector(Kappa(patchi));
}


Foam::tmp<Foam::surfaceScalarField> Foam::solidBlastThermo::q() const
{
    const fvMesh& mesh = this->T_.mesh();
    mesh.setFluxRequired(this->T_.name());

    return
      - (
            isotropic()
          ? fvm::laplacian(this->kappa(), this->T_)().flux()
          : fvm::laplacian(KappaLocal(), this->T_)().flux()
        )/mesh.magSf();
}


Foam::tmp<Foam::fvScalarMatrix> Foam::solidBlastThermo::divq
(
    volScalarField& e
) const
{
    return
      - (
            isotropic()
          ?   fvc::laplacian(this->kappa(), this->T_)
            + correction(fvm::laplacian(this->alpha(), e))
          :   fvc::laplacian(KappaLocal(), this->T_)
            + correction
              (
                  fvm::laplacian
                  (
                      KappaLocal()/this->Cv(),
                      e,
                      "laplacian(" + this->alpha().name() + ",e)"
                  )
              )
        );
}


// ************************************************************************* //
