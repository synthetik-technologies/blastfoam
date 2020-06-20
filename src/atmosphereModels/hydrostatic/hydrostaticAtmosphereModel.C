/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "hydrostaticAtmosphereModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace atmosphereModels
{
    defineTypeNameAndDebug(hydrostatic, 0);
    addToRunTimeSelectionTable(atmosphereModel, hydrostatic, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmosphereModels::hydrostatic::hydrostatic
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    atmosphereModel(mesh, dict),
    rhoRef_("rhoRef", dimDensity, dict_),
    pRef_("pRef", dimPressure, dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::atmosphereModels::hydrostatic::~hydrostatic()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmosphereModels::hydrostatic::createAtmosphere
(
    volScalarField& p,
    volScalarField& rho
) const
{
    rho = rhoRef_;
    p = pRef_ + (g_ & normal_)*rho*(h_ - groundElevation_);
}

// ************************************************************************* //
