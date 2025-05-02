/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "centrifugalBodyForce.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "boundBox.H"
#include "OStringStream.H"
#include "IStringStream.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(centrifugalBodyForce, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        centrifugalBodyForce,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::centrifugalBodyForce::centrifugalBodyForce
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, t, dict),
    angularVelocity_(inv(dimTime), dict.lookup("angularVelocity"))
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::centrifugalBodyForce::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    dict.readIfPresent("angularVelocity", angularVelocity_);

    return true;
}


bool Foam::functionObjects::centrifugalBodyForce::execute()
{
    return store
    (
        "centrifugalBodyForce",
        -(angularVelocity_ ^ (angularVelocity_ ^ mesh_.C()))
    );
}


bool Foam::functionObjects::centrifugalBodyForce::write()
{
    return writeObject("centrifugalBodyForce");
}

// ************************************************************************* //
