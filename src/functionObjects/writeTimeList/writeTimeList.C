/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies: | Added writeTimeList functionObject
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

#include "writeTimeList.H"
#include "SortableList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(writeTimeList, 0);
    addToRunTimeSelectionTable(functionObject, writeTimeList, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::writeTimeList::writeTimeList
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeTimes_(0),
    index_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::writeTimeList::~writeTimeList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::writeTimeList::read
(
    const dictionary& dict
)
{
    writeTimes_ = dict.lookup<scalarList>("times");
    writeTimes_.append(great);

    // Sort times from smallest to largest
    scalarList sortedTimes(writeTimes_);
    sort(sortedTimes);

    // Remove duplicate entries
    label j = 0;
    forAll(sortedTimes, i)
    {
        while (sortedTimes[i] == sortedTimes[i+1])
        {
            i++;
        }
        writeTimes_[j++] = sortedTimes[i];
    }
    writeTimes_.resize(j);

    // Get current location
    forAll(writeTimes_, ti)
    {
        if (writeTimes_[ti] > obr_.time().value())
        {
            index_ = ti;
            break;
        }
    }

    return true;
}


bool Foam::functionObjects::writeTimeList::execute()
{
    return true;
}


Foam::scalar Foam::functionObjects::writeTimeList::timeToNextWrite()
{
    return max(writeTimes_[index_] - obr_.time().value(), 0.0);
}


bool Foam::functionObjects::writeTimeList::write()
{
    if (mag(this->timeToNextWrite()) == 0)
    {
        Time& time(const_cast<Time&>(time_));
        time.writeNow();
        index_++;
    }
    return true;
}


// ************************************************************************* //
