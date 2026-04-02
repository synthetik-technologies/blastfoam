/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
     \\/     M anipulation  | Synthetik Applied Technologies
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

#include "timeIntegrationSystem.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeIntegrationSystem::timeIntegrationSystem(const word& name)
:
    timeIntegrationSystemBase(name),
    meshPtr_(nullptr),
    fvTimeInt_(nullptr),
    storeOld_(true)
{}

Foam::timeIntegrationSystem::timeIntegrationSystem
(
    const word& name,
    const fvMesh& mesh
)
:
    timeIntegrationSystemBase(name, mesh),
    meshPtr_(&mesh),
    fvTimeInt_
    (
        this->timeInt_.valid()
      ? dynamic_cast<const fvTimeIntegrator*>(this->timeInt_.ptr())
      : nullptr
    ),
    storeOld_(true)
{}


void Foam::timeIntegrationSystem::set(const fvMesh& mesh)
{
    timeIntegrationSystemBase::set(mesh);
    if (!meshPtr_.valid())
    {
        meshPtr_.reset(&mesh);
    }
    if (!fvTimeInt_.valid())
    {
        fvTimeInt_.reset
        (
            dynamic_cast<const fvTimeIntegrator*>(this->timeInt_.ptr())
        );
    }
}


void Foam::timeIntegrationSystem::set(const timeIntegrator& integrator)
{
    timeIntegrationSystemBase::set(integrator);
    if (!fvTimeInt_.valid())
    {
        fvTimeInt_.reset(dynamic_cast<const fvTimeIntegrator*>(&integrator));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::timeIntegrationSystem::~timeIntegrationSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
