/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Added writeTimeList functionObject
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
    writeTimes_(dict.lookup("times")),
    index_(0)
{
    writeTimes_.append(great);
    SortableList<scalar> sortedTimes(writeTimes_);
    sortedTimes.sort();
    writeTimes_ = sortedTimes;

    forAll(writeTimes_, ti)
    {
        if (writeTimes_[ti] > runTime.value())
        {
            index_ = ti;
            break;
        }
    }
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
    writeTimes_ = dict.lookupType<scalarList>("times");
    writeTimes_.append(great);
    SortableList<scalar> sortedTimes(writeTimes_);
    sortedTimes.sort();
    writeTimes_ = sortedTimes;

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
    return writeTimes_[index_] - obr_.time().value();
}


bool Foam::functionObjects::writeTimeList::write()
{
    if (mag(this->timeToNextWrite()) < small)
    {
        Time& time(const_cast<Time&>(time_));
        time.writeNow();
        index_++;
    }
    return true;
}


// ************************************************************************* //
