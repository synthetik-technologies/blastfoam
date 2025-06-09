/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2024-2025
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "hydrostaticPressureBoundaryTraction.H"
#include "uniformDimensionedFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace boundaryTractions
{
    defineTypeNameAndDebug(hydrostaticPressure, 0);
    addToRunTimeSelectionTable
    (
        boundaryTraction,
        hydrostaticPressure,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryTractions::hydrostaticPressure::hydrostaticPressure
(
    const word& name,
    const dictionary& dict,
    const feMesh1& femesh,
    const pointVectorField* DPtr
)
:
    boundaryTraction(name, dict, femesh, DPtr),
    pRef_(dict.lookup<scalar>("pRef")),
    hRef_(dict.lookupOrDefault<scalar>("hRef", 0.0)),
    rho_(dict.lookup<scalar>("rho")),
    g_
    (
        mesh_.foundObject<uniformDimensionedVectorField>("g")
      ? mesh_.lookupObject<uniformDimensionedVectorField>("g").value()
      : dict.lookup<vector>("g")
    ),
    magg_(mag(g_))
{}


Foam::boundaryTractions::hydrostaticPressure::hydrostaticPressure
(
    const hydrostaticPressure& hpbt
)
:
    boundaryTraction(hpbt),
    pRef_(hpbt.pRef_),
    hRef_(hpbt.hRef_),
    rho_(hpbt.rho_),
    g_(hpbt.g_),
    magg_(hpbt.magg_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryTractions::hydrostaticPressure::~hydrostaticPressure()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::boundaryTractions::hydrostaticPressure::pressure
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
    const UIndirectList<vector> nodes
    (
        femesh_.boundary()[patchID_].localNodes(),
        nodeLabels
    );
    scalar gh = hRef_*magg_;

    vector x(Zero);
    if (DPtr_)
    {
        const labelList& meshNodes =
            femesh_.boundary()[patchID_].meshNodes();
        const vectorField& D = *DPtr_;
        forAll(nodes, i)
        {
            gh += ((nodes[i] + D[meshNodes[nodeLabels[i]]]) & g_)*shape[i];
        }
    }
    else
    {
        forAll(nodes, i)
        {
            gh += (nodes[i] & g_)*shape[i];
        }
    }
    return pRef_ + rho_*gh;
}


Foam::vector Foam::boundaryTractions::hydrostaticPressure::traction
(
    const labelList& nodeLabels,
    const scalarList& shape,
    const vector& n
) const
{
   return Zero;
}


// ************************************************************************* //
