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
#include "UautoPtr.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Private Functinos * * * * * * * * * * * * * * //

template<class Type>
const Foam::dictionary& Foam::readComponent
(
    const dictionary& parentDict,
    const word& name,
    autoPtr<Modifier<Type>>& mod,
    List<Type>& values,
    const bool canRead
)
{
    UautoPtr<const entryTable> tablePtr;
    UautoPtr<const dictionary> dictPtr;
    if (parentDict.found(name + "File") && canRead)
    {
        tablePtr.set
        (
            &read2DTable
            (
                parentDict.lookup<fileName>(name + "File"),
                parentDict.lookupOrDefault<token>(name + "Delim", token::COMMA).pToken(),
                parentDict.lookupOrDefault<label>(name + "StartRow", 0),
                parentDict.lookupOrDefault<Switch>(name + "FlipTable", false)
            )
        );
    }
    else if (parentDict.found("file") && canRead)
    {
        tablePtr.set
        (
            &read2DTable
            (
                parentDict.lookup<fileName>("file"),
                parentDict.lookupOrDefault<token>("delim", token::COMMA).pToken(),
                parentDict.lookupOrDefault<label>("startRow", 0),
                parentDict.lookupOrDefault<Switch>("flipTable", false)
            )
        );
    }

    bool readFromTable = tablePtr.valid();

    label col = -1;
    label row = -1;
    label scale = 1;
    if (parentDict.found(name + "Coeffs"))
    {
        const dictionary& dict(parentDict.subDict(name + "Coeffs"));
        dictPtr.set(&dict);
        mod = Modifier<Type>::New
        (
            dict.lookupOrDefault<word>("mod", "none"),
            dict
        );
        mod->readReal(dict, "isReal");

        if (dict.found("scale"))
        {
            scale = dict.lookup<scalar>("scale");
        }

        if (!canRead)
        {}
        else if (dict.found(name))
        {
            values = dict.lookup<List<Type>>(name);
            readFromTable = false;
        }
        else if (dict.found("file"))
        {
            tablePtr.reset
            (
                &read2DTable
                (
                    dict.lookup<fileName>("file"),
                    dict.lookupOrDefault<token>("delim", token::COMMA).pToken(),
                    dict.lookupOrDefault<label>("startRow", 0),
                    dict.lookupOrDefault<Switch>("flipTable", false)
                )
            );

            if (tablePtr->m() == 1)
            {
                row = 0;
            }
            else if (tablePtr->n() == 1)
            {
                col = 0;
            }
            else
            {
                if (dict.found("col"))
                {
                    col = dict.lookup<label>("col");
                }
                else if (dict.found("row"))
                {
                    row = dict.lookup<label>("row");
                }
                else
                {
                    FatalIOErrorInFunction(dict)
                        << "Looking up a component of a 2D table requires either" << nl
                        << "a row (\"row\") or column (\"col\") to be specified" << endl
                        << abort(FatalIOError);
                }
            }

            readFromTable = true;
        }
        else if (dict.found("n"))
        {
            label ny = dict.lookup<label>("n");
            Type miny = dict.lookup<Type>("min");
            Type dy(Zero);
            if (!canRead)
            {}
            if (dict.found("delta"))
            {
                dy = dict.lookup<Type>("delta");
            }
            else if (dict.found("max"))
            {
                dy = (dict.lookup<Type>("max") - miny)/scalar(ny);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Either delta or max must be provided" <<endl
                    << abort(FatalIOError);
            }

            values.resize(ny);
            forAll(values, j)
            {
                values[j] = miny + dy*j;
            }
            readFromTable = false;
        }
        else if (!readFromTable)
        {
            FatalIOErrorInFunction(dict)
                << "Could not determine construction method of " << name << nl
                << "Pease provide a file to read from or (n, min, delta/max)" << endl
                << abort(FatalIOError);
        }

        if (readFromTable && col < 0 && row < 0)
        {
            if (dict.found("col"))
            {
                col = dict.lookup<label>("col");
            }
            else if (dict.found("row"))
            {
                row = dict.lookup<label>("row");
            }
        }
    }
    else if (parentDict.found("n" + name.capitalise()))
    {
        mod = Modifier<Type>::New
        (
            parentDict.lookupOrDefault<word>(name + "Mod", "none"),
            parentDict
        );
        mod->readReal(parentDict, name + "IsReal");

        if (canRead)
        {
            label ny = parentDict.lookup<label>("n" + name.capitalise());
            Type miny = parentDict.lookup<Type>("min" + name.capitalise());
            Type dy(Zero);
            if (parentDict.found("delta" + name.capitalise()))
            {
                dy = parentDict.lookup<Type>("delta" + name.capitalise());
            }
            else if (parentDict.found("max" + name.capitalise()))
            {
                dy =
                    (
                        parentDict.lookup<Type>("max" + name.capitalise())
                    - miny
                    )/scalar(ny);
            }
            else
            {
                FatalIOErrorInFunction(parentDict)
                    << "Either delta" << name.capitalise()
                    << " or max" << name.capitalise() << " must be provided" <<endl
                    << abort(FatalIOError);
            }

            values.resize(ny);
            forAll(values, j)
            {
                values[j] = miny + dy*j;
            }
            readFromTable = false;
        }
    }
    else if (parentDict.found(name) || readFromTable || !canRead)
    {
        if (parentDict.found(name))
        {
            values = parentDict.lookup<List<Type>>(name);
            readFromTable = false;
        }

        mod = Modifier<Type>::New
        (
            parentDict.lookupOrDefault<word>(name + "Mod", "none"),
            parentDict
        );
        mod->readReal(parentDict, name + "IsReal");
    }
    else
    {
        FatalIOErrorInFunction(parentDict)
            << "Neither the entry \"" << name << "\", "
            << "construction method "
            << "(n" << name.capitalise() << ", "
            << "min" << name.capitalise() << ", "
            << "max" << name.capitalise()
            << "/delta" << name.capitalise() << ") "
            << ", or the \""
            << name << "Coeffs\" subDictionary was found" << endl
            << abort(FatalIOError);
    }

    if (readFromTable)
    {
        const entryTable& table = tablePtr();
        if (col >= 0 || row >= 0)
        {}
        else if (table.m() == 1)
        {
            row = 0;
        }
        else if (table.n() == 1)
        {
            col = 0;
        }
        else
        {
            string colName(name + "Col");
            string rowName(name + "Row");
            if (parentDict.found(colName))
            {
                col = parentDict.lookup<label>(colName);
            }
            else if (parentDict.found(rowName))
            {
                row = parentDict.lookup<label>(rowName);
            }
            else
            {
                FatalIOErrorInFunction(parentDict)
                    << "Looking up a component of a 2D table requires either" << nl
                    << "a row (" << rowName << ") or column (" << colName << ")" << nl
                    << " to be specified" << endl
                    << abort(FatalIOError);
            }
        }

        if (col >= 0)
        {
            values = readColumn<Type>(table, col);
        }
        else if (row >= 0)
        {
            values = readRow<Type>(table, row);
        }
        else
        {
            FatalIOErrorInFunction(parentDict)
                << "Could not determine how to read table" << endl
                << abort(FatalIOError);
        }
    }
    if (canRead)
    {
        if (scale != 1.0)
        {
            values = scale*values;
        }
        if (!mod->isReal())
        {
            mod->Inv(values);
            mod->setReal();
        }
    }

    return dictPtr.valid() ? dictPtr() : parentDict;
}


template<class Type>
void Foam::read1DTable
(
    const fileName& file,
    const token::punctuationToken delim,
    List<Type>& values,
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
    const token::punctuationToken delim,
    List2D<Type>& data,
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
    DynamicList<List<Type>> tdata(10);
    while (is.good())
    {
        string line;
        is.getLine(line);
        removeComments(line);

        string lineEntry = line;
        lineEntry.replaceAll(delim, token::SPACE);
        lineEntry = '(' + lineEntry + ')';
        IStringStream iss(lineEntry);

        List<Type> lineVals(iss);

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
        if (data.m() != nx)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the x direction." << nl
                << "    Input: " << data.m() << nl
                << "    Read: " << nx << nl
                << abort(FatalError);
        }
        if (data.n() != ny)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the y direction." << nl
                << "    Input: " << data.n() << nl
                << "    Read: " << ny << nl
                << abort(FatalError);
        }
    }

    if (!flip)
    {
        data = tdata;
        return;
    }
    data.resize(nx, ny);
    for (label i = 0; i < data.m(); i++)
    {
        for (label j = 0; j < data.n(); j++)
        {
            data(i, j) = tdata[j][i];
        }
    }
}


template<class Type>
Foam::List<Type> Foam::readColumn
(
    const entryTable& entries,
    const label col
)
{
    if (!entries.size())
    {
        FatalErrorInFunction
            << "Trying to read from a table but no data was read." << endl
            << abort(FatalError);
    }
    if (col >= entries.n())
    {
        FatalErrorInFunction
            << "Only " << entries.n() << " columns were read, "
            << "but column " << col << " was requested" << endl
            << abort(FatalError);
    }

    List<Type> vals(entries.m());
    Type v;
    forAll(vals, i)
    {
        ITstream("column", entries(i, col))() >> v;
        vals[i] = v;
    }
    return vals;
}


template<class Type>
Foam::List<Type> Foam::readRow
(
    const entryTable& entries,
    const label row
)
{
    if (!entries.size())
    {
        FatalErrorInFunction
            << "Trying to read from a table but no data was read." << endl
            << abort(FatalError);
    }
    if (row >= entries.m())
    {
        FatalErrorInFunction
            << "Only " << entries.m() << " columns were read, "
            << " but row " << row << " was requested" << endl
            << abort(FatalError);
    }

    List<Type> vals(entries.n());
    Type v;
    forAll(vals, j)
    {
        ITstream("row", entries(row, j))() >> v;
        vals[j] = v;
    }
    return vals;
}


template<class Type>
void Foam::read3DTable
(
    const fileName& file,
    const token::punctuationToken delim,
    const token::punctuationToken rowDelim,
    List3D<Type>& data,
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

    List<List<List<Type>>> tdata;
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

        if (!strings.size())
        {
            continue;
        }
        else if (ny < 0)
        {
            ny = strings.size();
        }
        else if (ny != strings.size())
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Inconsistent size in the second direction" << endl
                << abort(FatalError);
        }

        List<List<Type>> yzvals(strings.size());
        forAll(strings, j)
        {
            string lineEntry = strings[j];
            lineEntry.replaceAll(delim, " ");
            lineEntry = '(' + lineEntry + ')';
            IStringStream iss(lineEntry);

            List<Type> zvals(iss);
            if (!yzvals.size())
            {
                continue;
            }
            else if (nz < 0)
            {
                nz = zvals.size();
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
        if (data.m() != nx)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the x direction." << nl
                << "    Input: " << data.m() << nl
                << "    Read: " << nx << nl
                << abort(FatalError);
        }
        if (data.n() != ny)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the y direction." << nl
                << "    Input: " << data.n() << nl
                << "    Read: " << ny << nl
                << abort(FatalError);
        }
        if (data.l() != nz)
        {
            FatalErrorInFunction
                << file << ":" << nl
                << "Size of input list is different that the size of the "
                << "read table in the z direction." << nl
                << "    Input: " << data.l() << nl
                << "    Read: " << nz << nl
                << abort(FatalError);
        }
    }
    else
    {
        data.resize(nx, ny, nz);
    }

    if (!flip)
    {
        data = tdata;
        return;
    }


    // Data needs to be correctly allocated before reading
    for (label i = 0; i < data.m(); i++)
    {
        for (label j = 0; j < data.n(); j++)
        {
            for (label k = 0; k < data.l(); k++)
            {
                data(i, j, k) = tdata[k][j][i];
            }
        }
    }
}


// ************************************************************************* //
