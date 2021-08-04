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

#include "phaseFluidBlastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(phaseFluidBlastThermo, 0);
    defineRunTimeSelectionTable(phaseFluidBlastThermo, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseFluidBlastThermo::phaseFluidBlastThermo
(
    const word& phaseName,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
:
    basicBlastThermo
    (
        phaseName,
        rho,
        e,
        T,
        dict,
        masterName
    ),
    mu_
    (
        lookupOrConstruct
        (
            rho.mesh(),
            IOobject::groupName("thermo:mu", phaseName),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            dimDynamicViscosity
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseFluidBlastThermo::~phaseFluidBlastThermo()
{}


// * * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseFluidBlastThermo> Foam::phaseFluidBlastThermo::New
(
    const word& phaseName,
    volScalarField& rho,
    volScalarField& e,
    volScalarField& T,
    const dictionary& dict,
    const word& masterName
)
{
    Info<<dictionaryConstructorTablePtr_->toc()<<endl;
    return basicBlastThermo::New<phaseFluidBlastThermo>
    (
        phaseName,
        rho,
        e,
        T,
        dict,
        masterName
    );
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::phaseFluidBlastThermo::mu() const
{
    return mu_;
}


Foam::tmp<Foam::scalarField>
Foam::phaseFluidBlastThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}


Foam::tmp<Foam::volScalarField> Foam::phaseFluidBlastThermo::nu() const
{
    return mu_/max(rho_, dimensionedScalar(dimDensity, 1e-10));
}


Foam::tmp<Foam::scalarField>
Foam::phaseFluidBlastThermo::nu(const label patchi) const
{
    return mu(patchi)/max(rho_.boundaryField()[patchi], 1e-10);
}


Foam::scalar Foam::phaseFluidBlastThermo::nui(const label celli) const
{
    return mu_[celli]/max(rho_[celli], 1e-10);
}

// ************************************************************************* //
