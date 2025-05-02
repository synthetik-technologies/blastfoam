/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2023
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

#include "commandEntry.H"
#include "OSspecific.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::commandEntry::commandEntry()
:
    command_("unknown"),
    options_(),
    nProcs_(1)
{}



Foam::commandEntry::commandEntry(Istream& is)
:
    command_("unknown"),
    options_(),
    nProcs_(1)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::commandEntry::~commandEntry()
{}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
Foam::Istream& Foam::operator>>(Istream& is, commandEntry& entry)
{
    string& str = entry;
    str.clear();
    entry.options_.clear();

    token t(is);
    if (t.isLabel())
    {
        entry.nProcs_ = t.labelToken();
        str = "mpirun -np " + Foam::name(entry.nProcs_) + " ";
    }
    else
    {
        entry.nProcs_ = 1;
        is.putBack(t);
    }

    string cmd(is);
    str += cmd;

    if (entry.nProcs_ > 1)
    {
        str += " -parallel";
    }

    IStringStream iss(cmd);
    entry.command_ = word(iss);

    iss >> t;
    while(iss.good())
    {
        if (!t.isPunctuation() || t.pToken() != token::SUBTRACT)
        {
            FatalIOErrorInFunction(iss)
                << "Expected " << token::SUBTRACT << " but found " << t
                << " while reading options" << endl
                << abort(FatalIOError);
        }
        {
            Tuple2<word, string> opt(word::null, string::null);
            iss >> t;
            opt.first() = t.wordToken();

            OStringStream oss;
            while
            (
                iss.read(t)
             && !t.isPunctuation()
             && t.pToken() != token::SUBTRACT
            )
            {
                oss << t << token::SPACE;
            }

            opt.second() = oss.str();
            if (opt.second().size())
            {
                opt.second().resize(opt.second().size()-1);
                entry.options_.append(opt);
            }
        }
        iss >> t;
    }
    return is;
}


// ************************************************************************* //
