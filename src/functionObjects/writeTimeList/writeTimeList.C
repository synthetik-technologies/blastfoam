/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Synthetik Applied Technologies: | Added writeTimeList functionObject
11-06-2025 Synthetik Applied Technologies: | Added spaced times
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

#include "IOmanip.H"
bool Foam::functionObjects::writeTimeList::read
(
    const dictionary& dict
)
{
    HashSet<scalar, Hash<scalar>> hashedTimes(writeTimes_);
    {
        ITstream is(dict.lookup("times"));
        token delim1, delim2, dtToken, nToken, tStartToken, tEndToken;
        token t(is);
        if (t.isLabel())
        {
            is >> t;
        }

        if (t.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "Incorrect first token, expected '(', found "
                << t.info()
                << exit(FatalIOError);
        }

        // Function to round time to the current time precision
        auto trimTime = [](const scalar t)
        {
            word tName(Time::timeName(t));
            OStringStream os;
            os << tName;
            IStringStream is(os.str());
            return readScalar(is);
        };


        while (is.good())
        {
            is >> t;
            if (t.isNumber())
            {
                tStartToken = t;

                is >> t;

                // Linear spacing with perscribed deltaT
                if (t.isPunctuation() && t.pToken() == token::COLON)
                {
                    delim1 = t;
                    is >> dtToken >> delim2 >> tEndToken;

                    if (!dtToken.isNumber())
                    {
                        FatalIOErrorInFunction(is)
                            << "Incorrect delta time token, expected number found "
                            << dtToken.info()
                            << exit(FatalIOError);
                    }
                    if (delim2.pToken() != token::COLON)
                    {
                        FatalIOErrorInFunction(is)
                            << "Incorrect second delimiter token, expected ':', found "
                            << delim2.info()
                            << exit(FatalIOError);
                    }
                    if (!tEndToken.isNumber())
                    {
                        FatalIOErrorInFunction(is)
                            << "Incorrect end time, number found "
                            << tEndToken.info()
                            << exit(FatalIOError);
                    }

                    scalar time = tStartToken.number();
                    const scalar dt = dtToken.number();
                    const scalar tEnd = tEndToken.number();
                    while (time <= tEnd)
                    {
                        hashedTimes.insert(trimTime(time));
                        time += dt;
                    }
                }
                else
                {
                    // Single entry
                    is.putBack(t);
                    hashedTimes.insert(trimTime(tStartToken.number()));
                }
            }
            // Function spacing
            else if (t.isWord())
            {
                word w(t.wordToken());
                bool isLog = w.find("log(") != string::npos;
                bool isExp = w.find("exp(") != string::npos;
                bool isLog10 = w.find("log10(") != string::npos;
                bool isPow10 = w.find("pow10(") != string::npos;
                bool isLinspace = w.find("linspace(") != string::npos;
                if (isLog)
                {
                    w.replace("log", "\0");
                }
                else if (isExp)
                {
                    w.replace("exp", "\0");
                }
                else if (isLog10)
                {
                    w.replace("log10", "\0");
                }
                else if (isPow10)
                {
                    w.replace("pow10", "\0");
                }
                else if (isLinspace)
                {
                    w.replace("linspace", "\0");
                }

                IStringStream iss(w);
                ITstream its("times", tokenList(iss));
                its >> tStartToken
                    >> delim1
                    >> nToken
                    >> delim2
                    >> tEndToken;
                if (!tStartToken.isNumber())
                {
                    FatalIOErrorInFunction(is)
                        << "Incorrect start time, number found "
                        << tStartToken.info()
                        << exit(FatalIOError);
                }
                if (delim1.pToken() != token::COLON)
                {
                    FatalIOErrorInFunction(is)
                        << "Incorrect first delimiter token, expected ':', found "
                        << delim1.info()
                        << exit(FatalIOError);
                }
                if (!nToken.isLabel())
                {
                    FatalIOErrorInFunction(is)
                        << "Incorrect token, expected label found "
                        << nToken.info()
                        << exit(FatalIOError);
                }
                if (delim2.pToken() != token::COLON)
                {
                    FatalIOErrorInFunction(is)
                        << "Incorrect second delimiter token, expected ':', found "
                        << delim2.info()
                        << exit(FatalIOError);
                }
                if (!tEndToken.isNumber())
                {
                    FatalIOErrorInFunction(is)
                        << "Incorrect end time, number found "
                        << tEndToken.info()
                        << exit(FatalIOError);
                }

                scalar tStart = tStartToken.number();
                const scalar nt = nToken.number();
                scalar tEnd = tEndToken.number();

                // Scale start and end times
                if (isLog)
                {
                    tStart = ::Foam::log(tStart);
                    tEnd = ::Foam::log(tEnd);
                }
                else if (isExp)
                {
                    tStart = ::Foam::exp(tStart);
                    tEnd = ::Foam::exp(tEnd);
                }
                else if (isLog10)
                {
                    tStart = log10(tStart);
                    tEnd = log10(tEnd);
                }
                else if (isPow10)
                {
                    tStart = pow(10.0, tStart);
                    tEnd = pow(10.0, tEnd);
                }

                scalar time = tStart;
                const scalar dt = (tEnd - tStart)/scalar(nt - 1);

                // Return to real time (inverse function)
                while (time <= tEnd)
                {
                    if (isLog)
                    {
                        hashedTimes.insert(trimTime(exp(time)));
                    }
                    else if (isExp)
                    {
                        hashedTimes.insert(trimTime(::Foam::log(time)));
                    }
                    else if (isLog10)
                    {
                        hashedTimes.insert(trimTime(pow(10, time)));
                    }
                    else if (isPow10)
                    {
                        hashedTimes.insert(trimTime(log10(time)));
                    }
                    else if (isLinspace)
                    {
                        hashedTimes.insert(trimTime(time));
                    }

                    time += dt;
                }
            }

            if (t.isPunctuation() && t.pToken() == token::END_LIST)
            {
                break;
            }
        }
    }

    // Return sorted times
    writeTimes_ = hashedTimes.sortedToc();
    writeTimes_.append(great);
    Info<<writeTimes_<<endl;
    std::exit(0);

    // Get current time index
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


Foam::scalar Foam::functionObjects::writeTimeList::timeToNextAction()
{
    return max(writeTimes_[index_] - obr_.time().value(), 0.0);
}


bool Foam::functionObjects::writeTimeList::write()
{
    if (mag(this->timeToNextAction()) < vSmall)
    {
        Time& time(const_cast<Time&>(time_));
        time.writeNow();
        index_++;
    }
    return true;
}


// ************************************************************************* //
