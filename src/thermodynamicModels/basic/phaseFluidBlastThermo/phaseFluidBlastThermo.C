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
    p_
    (
        rho.mesh().lookupObject<volScalarField>
        (
            IOobject::groupName("p", masterName)
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

const Foam::volScalarField& Foam::phaseFluidBlastThermo::p() const
{
    return p_;
}

// ************************************************************************* //
