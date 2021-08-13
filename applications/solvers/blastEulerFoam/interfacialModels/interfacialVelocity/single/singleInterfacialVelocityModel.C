/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "singleInterfacialVelocityModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialVelocityModels
{
    defineTypeNameAndDebug(single, 0);
    addToRunTimeSelectionTable(interfacialVelocityModel, single, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfacialVelocityModels::single::single
(
    const dictionary& dict,
    const phasePair& pair
)
:
    interfacialVelocityModel(dict, pair),
    phaseName_(dict.lookup("phase")),
    phase_
    (
        pair_.phase1().name() == phaseName_
      ? pair_.phase1()
      : pair_.phase2()
    ),
    U_(phase_.U()),
    phi_(phase_.phi())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interfacialVelocityModels::single::~single()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volVectorField>
Foam::interfacialVelocityModels::single::UI() const
{
    return U_;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::interfacialVelocityModels::single::phiI() const
{
    return phi_;
}

// ************************************************************************* //
