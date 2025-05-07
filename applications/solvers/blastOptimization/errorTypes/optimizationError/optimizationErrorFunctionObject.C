/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2025 Synthetik Applied Technologies: | Calculate optimizationError
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

#include "optimizationErrorFunctionObject.H"
#include "OFstream.H"
#include "OSspecific.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(optimizationError, 0);
    addToRunTimeSelectionTable(functionObject, optimizationError, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::optimizationError::optimizationError
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict)
{
    IOdictionary optimizationProperties
    (
        IOobject
        (
            "optimizationProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    errors_ = PtrList<errorType>
    (
        optimizationProperties.lookup("errors"),
        errorType::iNew(runTime, mesh_.name())
    );
    optimizationProperties.lookup("results") >> output_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::optimizationError::~optimizationError()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::optimizationError::read
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::functionObjects::optimizationError::execute()
{
    // Update errors
    Info<< "Current optimization status:" << incrIndent << endl;

    forAll(errors_,i)
    {
        errors_[i].update();
        Info<< indent
            << errors_[i].name() << " = " << errors_[i].value() << endl;
    }
    Info<< decrIndent << endl;

    return true;
}


bool Foam::functionObjects::optimizationError::write()
{
    mkDir(output_.path());

    OFstream os(output_);

    scalar totalError = 0.0;
    forAll(errors_, i)
    {
        const scalar errori = errors_[i].error();

        os << errors_[i].name() << ' ' << errors_[i].value() << ' ' << errori << nl;
        totalError += errors_[i].weight()*errori;
    }
    os << "totalError" << ' '<< totalError << ' ' << totalError << flush;
    return os.good();
}

// ************************************************************************* //
