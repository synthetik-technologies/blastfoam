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
    defineTypeNameAndDebug(centrifugalBodyForce, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        centrifugalBodyForce,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::centrifugalBodyForce::setBodyForce()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(polyMesh::defaultRegion);

    volVectorField& bodyForce =
        const_cast<volVectorField&>
        (
            mesh.lookupObject<volVectorField>("bodyForce")
        );

    bodyForce = -(angularVelocity_^(angularVelocity_^mesh.C()));

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::centrifugalBodyForce::centrifugalBodyForce
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    angularVelocity_(dict.lookup("angularVelocity"))
{
    Info << "Creating " << this->name() << " function object." << endl;

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::centrifugalBodyForce::start()
{
    return setBodyForce();
}


bool Foam::centrifugalBodyForce::execute()
{
    return setBodyForce();
}


bool Foam::centrifugalBodyForce::read(const dictionary& dict)
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return true;
}


bool Foam::centrifugalBodyForce::write()
{
    return setBodyForce();
}

// ************************************************************************* //
