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

#include "optimizationEntry.H"
#include "OSspecific.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::HashTable<Foam::dictionary> Foam::optimizationEntry::dictionaries;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::token* Foam::optimizationEntry::findTokenEntry
(
    const dictionary& dict,
    const fileName& path,
    label index
)
{
    index++;

    // Get components of the path
    List<wordRe> cmpts(path.components());
    forAll(cmpts, i)
    {
        cmpts[i].compile(wordRe::compOption::detect);
    }
    // Move up the path until no more dictionaries match
    const dictionary* dictPtr = &dict;
    label cmpti = 0;
    while (dictPtr->isDict(cmpts[cmpti]))
    {
        dictPtr = &dictPtr->subDict(cmpts[cmpti]);
        cmpti++;
        if (cmpti == cmpts.size()-1)
        {
            //If this reaches the end of the path, stop
            break;
        }
    }
    // Lookup the next component of the path
    const entry& e = dictPtr->lookupEntry(cmpts[cmpti++], true, true);
    if (cmpti == cmpts.size())
    {
        return const_cast<token*>(&(e.stream()[index]));
    }

    // Get the last part since that is the actual keyword and resize
    forAll(e.stream(), i)
    {
        // Only check strings
        if (e.stream()[i].isAnyString())
        {
            wordRe key
            (
                e.stream()[i].anyStringToken(),
                wordRe::compOption::detect
            );
            if (key.match(cmpts[cmpti]))
            {
                if (cmpti == cmpts.size()-1)
                {
                    return const_cast<token*>(&(e.stream()[index + i + 1]));
                }
                cmpts[cmpti++];
            }
        }
    }

    FatalIOErrorInFunction(dict)
        << "Could not find entry " << path << endl
        << abort(FatalIOError);
    return nullptr;
}


Foam::fileName Foam::optimizationEntry::findDictFromPath(fileName& path)
{
    label cmpti = 0;
    wordList cmpts(path.components());
    fileName file = cmpts[0];
    if
    (
        file == "constant"
        || file == "system"
    )
    {
        for (cmpti = 1; cmpti < cmpts.size(); cmpti++)
        {
            file = file / cmpts[cmpti];
            if (exists(file))
            {
                break;
            }
        }
    }
    else if (exists("constant" / file))
    {
        file = "constant" / file;
    }
    else if (exists("system" / file))
    {
        file = "system" / file;
    }

    if (!dictionaries.found(file))
    {
        fileName origFile(file + ".orig");

        if (exists(origFile) && exists(file, false))
        {
            rm(file);
        }
        // Copy file to *.orig to preserve original file
        else if (!exists(origFile, false))
        {
            cp(file, origFile);
        }

        IFstream ifs(file);
        if (!ifs.good())
        {
            FatalIOErrorInFunction(ifs)
                << "Could not open " << file << endl
                << abort(FatalIOError);
        }
        const bool disableFunctionEntries = entry::disableFunctionEntries;
        entry::disableFunctionEntries = true;
        dictionaries.insert(file, dictionary(ifs, true));
        entry::disableFunctionEntries = disableFunctionEntries;

    }

    path.clear();
    for (label cmptj = cmpti+1; cmptj < cmpts.size(); cmptj++)
    {
        path = path / cmpts[cmptj];
    }

    return file;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::optimizationEntry::optimizationEntry()
:
    file_("unknown"),
    path_("unknown"),
    index_(-1),
    t_(nullptr)
{}



Foam::optimizationEntry::optimizationEntry(Istream& is)
:
    file_("unknown"),
    path_("unknown"),
    index_(-1),
    t_(nullptr)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::optimizationEntry::~optimizationEntry()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::optimizationEntry::read(Istream& is)
{
    is >> path_;

    token t(is);
    if (t.isString())
    {
        file_ = path_;
        path_ = t.stringToken();
    }
    else
    {
        file_ = optimizationEntry::findDictFromPath(path_);
        is.putBack(t);
    }

    is >> t;
    if (t.isPunctuation() && t.pToken() == token::BEGIN_SQR)
    {
        is >> t;
        index_ = t.labelToken();
        is >> t;
        if (!t.isPunctuation() || t.pToken() != token::END_SQR)
        {
            FatalIOErrorInFunction(is)
                 << "Expected a '" << token::END_SQR
                << ", found " << t.info()
                << exit(FatalIOError);
        }
    }
    else
    {
        index_ = -1;
        is.putBack(t);
    }

    // Find the token and return the pointer
    t_ =
        optimizationEntry::findTokenEntry
        (
            optimizationEntry::dictionaries[file_],
            path_,
            index_
        );
}

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //
Foam::Istream& Foam::operator>>(Istream& is, optimizationEntry& entry)
{
    entry.read(is);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const optimizationEntry& entry)
{
    os  << entry.file_ << ": "
        << word(entry.path_);
    if (entry.index_ >= 0)
    {
        os  << " (" << entry.index_ << ")";
    }
    os  << endl;
    return os;
}



// ************************************************************************* //
