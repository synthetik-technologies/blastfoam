/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
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

#include "equationBase.H"
#include "OStringStream.H"

// * * * * * * * * * * Protected Static Member Functions * * * * * * * * * * //

Foam::string Foam::equationBase::mergeStrings(const List<string>& eqns) const
{
    if (!eqns.size())
    {
        return string::null;
    }
    OStringStream os;
    os << word(eqns[0]);
    for (label i = 1; i < eqns.size(); i++)
    {
        os << nl << word(eqns[i]);
    }
    return os.str();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equationBase::equationBase()
:
    equationBase(string::null)
{}


Foam::equationBase::equationBase(const string& eqnString)
:
    log_(false),
    logFile_(fileName::null),
    append_(false),
    logPtr_(nullptr),
    eqnString_(eqnString)
{}


Foam::equationBase::equationBase(const List<string>& eqnStrings)
:
    equationBase(mergeStrings(eqnStrings))
{}


Foam::equationBase::equationBase(const string& eqnString, const dictionary& dict)
:
    log_(dict.lookupOrDefault("log", false)),
    logFile_
    (
        dict.lookupOrDefault<fileName>
        (
            "logFile",
            "${FOAM_CASE}.evals"
        )
    ),
    append_(dict.lookupOrDefault("append", false)),
    logPtr_(nullptr),
    eqnString_(dict.lookupOrDefault<string>("eqnString", eqnString))
{
    const fileName origLogFile(logFile_);
    logFile_.expand();
}


Foam::equationBase::equationBase
(
    const List<string>& eqnStrings,
    const dictionary& dict
)
:
    equationBase(mergeStrings(eqnStrings), dict)
{}


Foam::equationBase::equationBase(const dictionary& dict)
:
    equationBase(word::null, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equationBase::~equationBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::equationBase::obr() const
{
    if (!obrPtr_.valid())
    {
        FatalErrorInFunction
            << "Trying to access equation objectRegistry, but it has not been set"
            << endl
            << abort(FatalError);
    }
    return obrPtr_();
}


const Foam::fileName& Foam::equationBase::logFileName() const
{
    return logFile_;
}


void Foam::equationBase::setLog(const fileName& logFile, const bool log)
{
    log_ = log;
    if (logPtr_.valid() && logPtr_->name() != logFile)
    {
        logPtr_.clear();
    }
    logFile_ = logFile;
}


void Foam::equationBase::setLog(const bool log)
{
    log_ = log;
}


Foam::OFstream& Foam::equationBase::logStream() const
{
    if (!logPtr_.valid())
    {
        logPtr_.set
        (
            new OFstream
            (
                logFile_,
                OFstream::ASCII,
                OFstream::currentVersion,
                OFstream::UNCOMPRESSED,
                append_
            )
        );
    }
    return logPtr_();
}

void Foam::equationBase::read(const dictionary& dict)
{
    dict.readIfPresent("log", log_);
    dict.readIfPresent("eqnString", eqnString_);

    if (log_)
    {
        if (dict.found("logFile") || logFile_.empty())
        {
            logFile_ = dict.lookup<fileName>("logFile");
            logFile_.expand();
        }
        dict.readIfPresent("append", append_);
    }
}

// ************************************************************************* //
