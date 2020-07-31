/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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

#include "afterburnModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(afterburnModel, 0);
    defineRunTimeSelectionTable(afterburnModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::afterburnModel::afterburnModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    integrationSystem
    (
        IOobject::groupName("afterburnModel", phaseName),
        mesh
    ),
    mesh_(mesh),
    dict_(dict),
    time_(mesh.time()),
    dt_(mesh.time().deltaT())
{
    this->lookupAndInitialize();
    times_.resize(this->oldIs_.size() + 1);
    forAll(times_, ti)
    {
        times_[ti].dimensions().reset(dimTime);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::afterburnModel::~afterburnModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::afterburnModel::solve
(
    const label stepi,
    const scalarList& ai,
    const scalarList& bi
)
{
    dt_ = mesh_.time().deltaT();
    dimensionedScalar t0(dimTime, 0.0);
    if (stepi == 1)
    {
        times_[0] = mesh_.time();
    }
    scalar f = 0;
    for (label i = 0; i < stepi; i++)
    {
        t0 += ai[i]*times_[i];
        f += bi[i];
    }
    dt_ *= f;
    time_ = t0 + dt_;
    times_[stepi] = time_;
}


// ************************************************************************* //
