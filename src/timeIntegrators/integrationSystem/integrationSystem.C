/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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

#include "integrationSystem.H"
#include "timeIntegrator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::integrationSystem::integrationSystem
(
    const word& name,
    const fvMesh& mesh
)
:
    mesh_(mesh),
    name_(name)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::integrationSystem::~integrationSystem()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::integrationSystem::lookupAndInitialize(const word& name)
{
    if (mesh_.foundObject<timeIntegrator>(name))
    {
        mesh_.lookupObjectRef<timeIntegrator>(name).setODEFields(*this);
    }
}


void Foam::integrationSystem::setODEFields
(
    const label nSteps,
    const boolList& storeFields,
    const boolList& storeDeltas
)
{
    oldIs_.resize(nSteps);
    deltaIs_.resize(nSteps);
    label fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeFields[i])
        {
            oldIs_[i] = fi;
            fi++;
        }
        else
        {
            oldIs_[i] = -1;
        }
    }
    nOld_ = fi;

    fi = 0;
    for (label i = 0; i < nSteps; i++)
    {
        if (storeDeltas[i])
        {
            deltaIs_[i] = fi;
            fi++;
        }
        else
        {
            deltaIs_[i] = -1;
        }
    }
    nDelta_ = fi;

    this->clearODEFields();
}


Foam::scalar Foam::integrationSystem::f
(
    const label stepi,
    const scalarList& b
) const
{
    scalar ft = b[stepi - 1];
    for (label j = 0; j < stepi - 1; j++)
    {
        label i = deltaIs_[j];
        if (i != -1)
        {
            ft += b[i];
        }
    }
    return ft;
}


bool Foam::integrationSystem::writeData(Ostream& os) const
{
    return os.good();
}

// ************************************************************************* //
