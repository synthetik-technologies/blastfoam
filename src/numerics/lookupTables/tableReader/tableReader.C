/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

void Foam::removeComments(string& line)
{
    label commentI = line.find('#');
    if (commentI >= 0)
    {
        line = line.substr(0, commentI);
    }
}


Foam::List<Foam::List<Foam::string>> Foam::read2DTable
(
    const fileName& file,
    const string& delim,
    const label startLine,
    const bool flip
)
{
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
    label lineI = 0;

    word line;
    DynamicList<List<string>> tentries;
    while (is.good())
    {
        is.getLine(line);
        if (lineI++ < startLine)
        {
            continue;
        }
        removeComments(line);

        DynamicList<word> lineVals(line.size());
        label stringi = 0;
        for
        (
            string::const_iterator iter = line.begin();
            iter != line.end();
            ++iter
        )
        {
            if (*iter != delim[0])
            {
                lineVals(stringi) = lineVals(stringi) + *iter;
            }
            else
            {
                stringi++;
            }
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
                << line
                << abort(FatalError);
        }
        tentries.append(List<string>(lineVals.size()));
        forAll(lineVals, i)
        {
            tentries[nx][i] = lineVals[i];
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


    if (!f)
    {
        return move(tentries);
    }
    List<List<string>> entries(nx, List<string>(ny));
    forAll(entries, i)
    {
        forAll(entries[i], j)
        {
            entries[i][j] = tentries[j][i];
        }
    }
    return entries;
}

// ************************************************************************* //
