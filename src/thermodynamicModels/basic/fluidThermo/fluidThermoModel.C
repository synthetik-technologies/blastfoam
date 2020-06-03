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

#include "fluidThermoModel.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(fluidThermoModel, 0);
    defineRunTimeSelectionTable(fluidThermoModel, basic);
    defineRunTimeSelectionTable(fluidThermoModel, detonating);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidThermoModel::fluidThermoModel
(
    const word& phaseName,
    volScalarField& p,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const bool master
)
:
    basicThermoModel
    (
        phaseName,
        p,
        rho,
        e,
        T,
        dict,
        master
    ),
    mu_
    (
        IOobject
        (
            IOobject::groupName("thermo:mu", phaseName),
            p_.mesh().time().timeName(),
            p_.mesh()
        ),
        p_.mesh(),
        dimensionedScalar(dimDynamicViscosity, 0.0),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0),
    viscous_(true)
{}


Foam::fluidThermoModel::fluidThermoModel
(
    const word& phaseName,
    const fvMesh& mesh,
    const dictionary& dict,
    const bool master
)
:
    basicThermoModel
    (
        phaseName,
        mesh,
        dict,
        master
    ),
    mu_
    (
        IOobject
        (
            IOobject::groupName("thermo:mu", phaseName),
            mesh.time().timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimDynamicViscosity, 0.0),
        wordList(p_.boundaryField().types().size(), "zeroGradient")
    ),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0),
    viscous_(true)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidThermoModel::~fluidThermoModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fluidThermoModel::initialize()
{
    if (!master_)
    {
        return;
    }
    if (max(e_).value() <= 0.0)
    {
        volScalarField e(calce());
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
        eBoundaryCorrection();
    }
    this->correct();

    if (max(mu_).value() < small)
    {
        viscous_ = false;
    }
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::fluidThermoModel::nu() const
{
    return mu_/rho_;
}


Foam::tmp<Foam::scalarField>
Foam::fluidThermoModel::nu(const label patchi) const
{
    return mu(patchi)/rho_.boundaryField()[patchi];
}


Foam::scalar Foam::fluidThermoModel::nui(const label celli) const
{
    return mu_[celli]/rho_[celli];
}

// ************************************************************************* //
