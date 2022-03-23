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

template<class Type>
bool Foam::readComponent
(
    const dictionary& parentDict,
    const word& name,
    word& modType,
    Field<Type>& values,
    const List<List<string>>& table
)
{
    Switch isReal = true;
    bool readFromTable = table.size();
    label col = -1;
    label scale = 1.0;
    if (parentDict.found(name + "Coeffs"))
    {
        const dictionary& dict(parentDict.subDict(name + "Coeffs"));
        modType = dict.lookupOrDefault<word>("mod", "none");

        if (dict.found("scale"))
        {
            scale = dict.lookup<scalar>("scale");
        }

        if (readFromTable)
        {
            if (dict.found("col"))
            {
                col = dict.lookup<label>("col");
            }
        }
        else if (dict.found(name))
        {
            values = dict.lookup<Field<Type>>(name);
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
        }
        else if (dict.found("file"))
        {
            fileName file(dict.lookup("file"));
            read1DTable
            (
                file,
                dict.lookupOrDefault<string>("delim", ","),
                values
            );
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
        }
        else
        {
            label ny = dict.lookup<label>("n");
            Type dy = dict.lookup<Type>("delta");
            Type miny = dict.lookup<Type>("min");

            values.resize(ny);
            forAll(values, j)
            {
                values[j] = miny + dy*j;
            }
            isReal = dict.lookupOrDefault<Switch>("isReal", true);
        }
    }
    else if (parentDict.found(name) || readFromTable)
    {
        if (readFromTable)
        {}
        else
        {
            values = parentDict.lookup<Field<Type>>(name);
            isReal =
                parentDict.lookupOrDefault<Switch>
                (
                    name + "isReal",
                    true
                );
        }

        modType = parentDict.lookupOrDefault<word>(name + "Mod", "none");
    }
    else
    {
        FatalErrorInFunction
            << "Neither the entry \"" << name << "\" or the \""
            << name << "Coeffs\" subDictionary was found" << endl
            << abort(FatalError);
    }

    if (readFromTable)
    {
        values = readColumn<Type>
        (
            table,
            col < 0 ? parentDict.lookup<label>(name + "Col") : col
        );
        isReal =
            parentDict.lookupOrDefault<Switch>
            (
                name + "isReal",
                true
            );
    }

    if (scale != 1.0)
    {
        values = scale*values;
    }
    return isReal;
}


template<class Type>
void Foam::read1DTable
(
    const fileName& file,
    const string& delim,
    Field<Type>& values,
    const bool determineSize
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

    List<Type> tvals;
    string line;
    while (is.good())
    {
        is.getLine(line);
        removeComments(line);
        line.replaceAll(delim, " ");
        line = '(' + line + ')';

        IStringStream isLine(line);
        List<Type> lineVals(isLine);
        if (lineVals.size())
        {
            tvals.append(lineVals);
        }
    }
    if (!determineSize && values.size() != tvals.size())
    {
        FatalErrorInFunction
            << file << ":" << nl
            << "Size of input list is different that the size of the read "
            << "list." << nl
            << "    Input: " << values.size() << nl
            << "    Read: " << tvals.size() << nl
            << abort(FatalError);
    }

    values = tvals;
}


template<class Type>
void Foam::read2DTable
(
    const fileName& file,
    const string& delim,
    Field<Field<Type>>& data,
    const bool flip,
    const bool determineSize
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
    Field<Field<Type>> tdata;
    while (is.good())
    {
        string line;
        is.getLine(line);
        removeComments(line);

        string lineEntry = line;
        lineEntry.replaceAll(delim, " ");
        lineEntry = '(' + lineEntry + ')';
        IStringStream iss(lineEntry);

        Field<Type> lineVals(iss);

        if (ny < 0)
        {
            ny = lineVals.size();
        }
        else if (!lineVals.size())
        {
            break;
        }
        else if (lineVals.size() != ny)
        {
            FatalErrorInFunction
                << "Incompatible table rows" << endl
                << line
                << abort(FatalError);
        }

        tdata.append(lineVals);
        nx++;
    }

    if (flip)
    {
        label t = nx;
        nx = ny;
        ny = t;
    }

    if (!determineSize)
    {
        if (data.size() != nx)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the x direction." << nl
                << "    Input: " << data.size() << nl
                << "    Read: " << nx << nl
                << abort(FatalError);
        }
        if (data[0].size() != ny)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the y direction." << nl
                << "    Input: " << data[0].size() << nl
                << "    Read: " << ny << nl
                << abort(FatalError);
        }
    }

    if (!flip)
    {
        data.transfer(tdata);
        return;
    }
    forAll(data, i)
    {
        forAll(data[i], j)
        {
            data[i][j] = tdata[j][i];
        }
    }
}


template<class Type>
Foam::List<Type> Foam::readColumn
(
    const List<List<string>>& entries,
    const label col
)
{
    if (!entries.size())
    {
        FatalErrorInFunction
            << "Trying to read from a table but no data was read." << endl
            << abort(FatalError);
    }
    if (col > entries.first().size())
    {
        FatalErrorInFunction
            << "Only " << entries.first().size() << " columns were read, " << nl
            << " but column " << col << " was requested" << endl
            << abort(FatalError);
    }

    List<Type> vals(entries.size());
    Type v;
    forAll(entries, i)
    {
        IStringStream(entries[i][col])() >> v;
        vals[i] = v;
    }
    return vals;
}


template<class Type>
void Foam::read3DTable
(
    const fileName& file,
    const string& delim,
    const string& rowDelim,
    Field<Field<Field<Type>>>& data,
    const bool flip,
    const bool determineSize
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

    Field<Field<Field<Type>>> tdata;
    label nx = 0;
    label ny = -1;
    label nz = -1;
    while (is.good())
    {
        string line;
        is.getLine(line);
        removeComments(line);

        label endSplit = line.find(rowDelim);
        stringList strings;
        while (line.size())
        {
            strings.append(line(0, endSplit));
            line = line(endSplit+1, line.size());
            endSplit = line.find(rowDelim);
        }

        if (ny < 0)
        {
            ny = strings.size();
        }
        else if (!strings.size())
        {
            continue;
        }
        else if (ny != strings.size())
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Inconsistent size in the second direction" << endl
                << abort(FatalError);
        }

        Field<Field<Type>> yzvals(strings.size());
        forAll(strings, j)
        {
            string lineEntry = strings[j];
            lineEntry.replaceAll(delim, " ");
            lineEntry = '(' + lineEntry + ')';
            IStringStream iss(lineEntry);

            Field<Type> zvals(iss);
            if (nz < 0)
            {
                nz = zvals.size();
            }
            else if (!yzvals.size())
            {
                continue;
            }
            else if (nz != zvals.size())
            {
                FatalErrorInFunction
                    << file << ":" << nl
                    << "Inconsistent size in the third direction" << endl
                    << abort(FatalError);
            }


            yzvals[j] = zvals;
        }
        tdata.append(yzvals);
        nx++;
    }

    if (flip)
    {
        label t = nx;
        nx = nz;
        nz = t;
    }

    if (!determineSize)
    {
        if (data.size() != nx)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the x direction." << nl
                << "    Input: " << data.size() << nl
                << "    Read: " << nx << nl
                << abort(FatalError);
        }
        if (data[0].size() != ny)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the y direction." << nl
                << "    Input: " << data[0].size() << nl
                << "    Read: " << ny << nl
                << abort(FatalError);
        }
        if (data[0][0].size() != nz)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the z direction." << nl
                << "    Input: " << data[0][0].size() << nl
                << "    Read: " << nz << nl
                << abort(FatalError);
        }
    }
    else
    {
        data.resize(nx);
        forAll(data, i)
        {
            data[i].resize(ny);
            forAll(data[i], j)
            {
                data[i][j].resize(nz);
            }
        }
    }

    if (!flip)
    {
        data.transfer(tdata);
        return;
    }


    // Data needs to be correctly allocated before reading
    forAll(data, i)
    {
        forAll(data[i], j)
        {
            forAll(data[i][j], k)
            {
                data[i][j][k] = tdata[k][j][i];
            }
        }
    }
}


// ************************************************************************* //
