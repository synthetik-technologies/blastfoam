/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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

#include "tableReader.H"


Foam::HashTable<Foam::entryTable> Foam::readTables;

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

void Foam::removeComments(string& line)
{
    label commentI = line.find('#');
    if (commentI >= 0)
    {
        line = line.substr(0, commentI);
    }
}


Foam::token::punctuationToken Foam::readDelim
(
    const dictionary& dict,
    const word& name,
    const token::punctuationToken delim
)
{
    if (!dict.found(name))
    {
        return delim;
    }
    ITstream is = dict.lookup(name);
    token t(is);
    if (!t.isString() || t.stringToken().size() != 1)
    {
        FatalIOErrorInFunction(is)
            << "Expected single quoted character but found " << t << endl
            << abort(FatalIOError);
    }
    return token::punctuationToken(t.stringToken()[0]);
}


const Foam::entryTable& Foam::read2DTable
(
    const fileName& file,
    const token::punctuationToken delim,
    const label startLine,
    const bool flip
)
{
    if (readTables.found(file))
    {
        return readTables[file];
    }

    fileName fNameExpanded(file);
    fNameExpanded.expand();

    // Open a stream and check it
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fNameExpanded));
    ISstream& is = isPtr();
    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open file" << file << nl
            << exit(FatalIOError);
    }

    DynamicList<Tuple2<scalar, scalar>> values;

    label ny = -1;
    label nx = 0;

    DynamicList<List<List<token>>> tentries;
    token t(is);

    while (is.good())
    {
        if (t.lineNumber() <= startLine)
        {
            is >> t;
            continue;
        }
        label lineNo = t.lineNumber();
        if
        (
            (t.isPunctuation() && t.pToken() == token::HASH)
         || t.isFunctionName()
        )
        {
            do
            {
                is >> t;
            } while
            (
                (t.isPunctuation() && t.pToken() != token::NL)
             && t.good()
             && t.lineNumber() == lineNo
            );
        }

        DynamicList<DynamicList<token>> lineVals;
        label cmpti = 0;
        lineNo = t.lineNumber();
        while (is.good())
        {
            bool add = true;
            if (!t.good())
            {
                break;
            }
            else if (t.lineNumber() != lineNo)
            {
                lineNo = t.lineNumber();
                break;
            }

            if (delim == token::SPACE)
            {
                if (lineVals(cmpti).size())
                {
                    cmpti++;
                }
            }
            else if (t.isPunctuation())
            {
                if (t.pToken() == token::NL)
                {
                    break;
                }
                if (t.pToken() == delim)
                {
                    if (lineVals(cmpti).size())
                    {
                        cmpti++;
                    }
                    add = false;
                }
            }
            if (add && t.good())
            {
                lineVals(cmpti).append(t);
            }
            is>> t;
        }

        if (!lineVals.size())
        {
            continue;
        }
        else if (ny < 0)
        {
            ny = lineVals.size();
        }
        else if (lineVals.size() != ny)
        {
            FatalErrorInFunction
                << "Incompatible table rows" << endl
                << abort(FatalError);
        }

        tentries.append(List<List<token>>(lineVals.size()));
        forAll(lineVals, i)
        {
            tentries[nx][i].transfer(lineVals[i]);
        }
        nx++;
    }

    // If only one row is provided, assume this is the data
    bool f = flip;
    if (flip || nx == 1)
    {
        f = true;
        label t = nx;
        nx = ny;
        ny = t;
    }

    entryTable& entries = readTables(file);
    entries.setSize(nx, ny);
    if (!f)
    {
        forAll(tentries, i)
        {
            forAll(tentries[i], j)
            {
                entries(i, j).transfer(tentries[i][j]);
            }
        }
    }
    else
    {
        forAll(tentries, j)
        {
            forAll(tentries[j], i)
            {
                entries(i, j).transfer(tentries[j][i]);
            }
        }
    }

    return entries;
}

// ************************************************************************* //
