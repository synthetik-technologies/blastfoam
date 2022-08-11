/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "massToCell.H"
#include "polyMesh.H"
#include "calcAngleFraction.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::massToCell::checkMass
(
    const labelHashSet& set,
    const polyMesh& mesh
) const
{
    if (mass_ < small)
    {
        return;
    }
    scalar sumV = 0.0;
    forAllConstIter(labelHashSet, set, iter)
    {
        sumV += mesh.cellVolumes()[iter.key()];
    }
    reduce(sumV, sumOp<scalar>());
    scalar actualMass = sumV*rho_/calcAngleFraction(mesh);
    if (mag(actualMass - mass_)/mass_ > 0.1)
    {
        WarningInFunction
            << "Requested mass is " << mass_
            << " but set mass is " << actualMass
            << ", " << mag(actualMass - mass_)/mass_*100.0 << "% different"
            << nl << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::massToCell::massToCell(const dictionary& dict)
:
    mustRead_(dict.dictName() != "backup"),
    rho_(1.0),
    mass_(0.0),
    centre_(dict.lookup<vector>("centre")),
    volume_(0.0),
    scale_(dict.lookupOrDefault("scale", 1.0)),
    read_(false)
{
    if (mustRead_)
    {
        rho_ = dict.lookup<scalar>("rho");
        mass_ = dict.lookup<scalar>("mass");
    }
    else
    {
        rho_ = dict.lookupOrDefault<scalar>("rho", 1.0);
        mass_ = dict.lookupOrDefault<scalar>("mass", 0.0);
    }
    volume_ = mass_/rho_*scale_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::massToCell::~massToCell()
{}


// ************************************************************************* //
