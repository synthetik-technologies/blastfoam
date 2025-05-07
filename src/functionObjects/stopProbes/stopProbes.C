/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025
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

#include "stopProbes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stopProbes, 0);
    addToRunTimeSelectionTable(functionObject, stopProbes, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stopProbes::stopProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldValues_(),
    locations_(),
    stopControl_(Time::stopAtControl::nextWrite),
    stop_(false)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::stopProbes::~stopProbes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stopProbes::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.lookup("locations") >> locations_;

    dict.lookup("fields") >> fieldValues_;

    if (dict.found("stopAt"))
    {
        stopControl_ = Time::stopAtControlNames.read(dict.lookup("stopAt"));
    }

    return true;
}


bool Foam::functionObjects::stopProbes::execute()
{
    if (!stop_)
    {
        mesh_.tetBasePtIs();
        labelList cells(locations_.size(), -1);
        forAll(locations_, i)
        {
            const point& pt = locations_[i];
            cells[i] = mesh_.findCell(pt);
        }

        forAllConstIter(HashTable<scalar>, fieldValues_, iter)
        {
            const volScalarField& fld =
                mesh_.lookupObject<volScalarField>(iter.key());
            const scalar val = iter();

            forAll(cells, i)
            {
                if (cells[i] >= 0 && fld[cells[i]] > val)
                {
                    stop_ = true;
                    break;
                }
            }
        }
        reduce(stop_, orOp<bool>());
    }

    if (stop_)
    {
        const_cast<Time&>(mesh_.time()).stopAt(stopControl_);
    }
    return true;
}

// ************************************************************************* //
