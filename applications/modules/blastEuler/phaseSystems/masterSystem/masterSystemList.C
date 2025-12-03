/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
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

#include "masterSystemList.H"
#include "masterSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(masterSystemList, 0);
}

Foam::masterSystemList& Foam::masterSystemList::New
(
    const objectRegistry& obr
)
{
    if (!obr.foundObject<masterSystemList>(masterSystemList::typeName))
    {
        masterSystemList* systemPtr
        (
            new masterSystemList(obr)
        );

        // Transfer ownership of this object to the objectRegistry
        systemPtr->store(systemPtr);
    }

    return obr.lookupObjectRef<masterSystemList>(masterSystemList::typeName);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::masterSystemList::masterSystemList
(
    const objectRegistry& obr
)
:
    regIOobject
    (
        IOobject
        (
            typeName,
            obr.time().constant(),
            obr
        )
    ),
    UPtrList<masterSystem>(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::masterSystemList::~masterSystemList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::masterSystemList::addSystem
(
    masterSystem& system
)
{
    const label i = this->size();
    this->resize(i + 1);
    this->set(i, &system);
}


void Foam::masterSystemList::initialize()
{
    forAll(*this, i)
    {
        (*this)[i].initialize();
    }
}


void Foam::masterSystemList::correctAlpha()
{
    forAll(*this, i)
    {
        (*this)[i].correctAlpha();
    }
}


void Foam::masterSystemList::update()
{
    forAll(*this, i)
    {
        (*this)[i].update();
    }
}


void Foam::masterSystemList::solve()
{
    forAll(*this, i)
    {
        (*this)[i].solve();
    }
}


void Foam::masterSystemList::postExplicit()
{
    forAll(*this, i)
    {
        (*this)[i].postExplicit();
    }
}


void Foam::masterSystemList::postImplicit()
{
    forAll(*this, i)
    {
        (*this)[i].postImplicit();
    }
}


void Foam::masterSystemList::storeExplicit()
{
    forAll(*this, i)
    {
        (*this)[i].storeExplicit();
    }
}


void Foam::masterSystemList::postUpdate()
{
    forAll(*this, i)
    {
        (*this)[i].postUpdate();
    }
}


void Foam::masterSystemList::clear()
{
    forAll(*this, i)
    {
        (*this)[i].clear();
    }
}

// ************************************************************************* //
