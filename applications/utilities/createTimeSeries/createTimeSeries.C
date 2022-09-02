/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2022
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

Description
    Utility to create a time series file for easier viewing of post processing
    data

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "PtrListDictionary.H"
#include "OFstream.H"
#include "argList.H"

using namespace Foam;

// Valid *.series files
const HashSet<word> validExts({"vtk", "stl"});

// Find relative path from source to destination
// source: a/b/c, destination: a/b/d/e, relPath: ../d/e
fileName relativePath(const fileName& source, const fileName& destination)
{
    // src:
    wordList srcCmpts(source.components());
    wordList destCmpts(destination.components());

    label starti = 0;
    for (label i = 0; i < srcCmpts.size(); i++)
    {
        if (srcCmpts[i] != destCmpts[i])
        {
            break;
        }
        starti++;
    }

    fileName relPath;
    if (starti != srcCmpts.size())
    {
        relPath = "..";
        for (label i = starti+1; i < srcCmpts.size(); i++)
        {
            relPath = relPath / "..";
        }
        relPath = relPath / destCmpts[starti++];
    }
    else
    {
        relPath = destCmpts[starti++];
    }
    for (label i = starti; i < destCmpts.size(); i++)
    {
        relPath = relPath / destCmpts[i];
    }
    return relPath;
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Creates a *.series file using the given path.\n"
        "Valid extensions are stl and vtk\n"
    );

    argList::validArgs.append("path");
    argList::addBoolOption
    (
        "absolute",
        "Create the time series using absolute file paths."
    );
    argList::addOption
    (
        "outputDir",
        "Directory to write outputs to"
    );

    #include "setRootCase.H"

    // Get the path to read files from
    fileName path(args.argRead<fileName>(1));
    fileName outputDir = args.optionLookupOrDefault("outputDir", path);
    outputDir.clean();
    mkDir(outputDir);

    Info<< "Reading valid time directories from " << path << endl;

    bool absolutePath = args.optionFound("absolute");
    if (absolutePath)
    {
        Info<< "Using absolute paths" << endl;
    }

    // Collect times
    fileNameList timeDirs
    (
        readDir
        (
            path,
            fileType::directory
        )
    );

    // Make sure that files are actually times
    DynamicList<scalar> times;
    forAll(timeDirs, ti)
    {
        token t((IStringStream(timeDirs[ti]))());
        if (t.isNumber())
        {
            times.append(t.number());
        }
    }

    // Sort times
    sort(times);
    timeDirs.setSize(times.size());
    forAll(times, ti)
    {
        timeDirs[ti] = Time::timeName(times[ti]);
    }

    //- List of *.series files
    PtrListDictionary<OFstream> seriesOutputs(0);

    Info<< endl;
    forAll(timeDirs, ti)
    {
        //- Check directory for valid files and create series outputs
        fileNameList files
        (
            readDir
            (
                path/timeDirs[ti],
                fileType::file
            )
        );

        forAll(files, i)
        {
            const word ext = files[i].ext();
            if (validExts.found(ext))
            {
                const fileName& file(files[i]);

                if (!seriesOutputs.found(file))
                {
                    Info<< "Making time series file for " << file << endl;

                    const fileName outputFile(outputDir/file + ".series");
                    Info<< "Writing " << outputFile << nl << endl;

                    seriesOutputs.resize(seriesOutputs.size() + 1);
                    seriesOutputs.set
                    (
                        seriesOutputs.size() - 1,
                        file,
                        new OFstream(outputFile)
                    );

                    // Header
                    seriesOutputs[file]
                        << '{' << nl << incrIndent
                        << indent << "\"file-series-version\" : \"1.0\"," << nl
                        << indent << "\"files\" : [" << incrIndent << nl;
                }
            }
        }
    }

    // Hash set to make sure comma is only written after first time is added
    HashSet<word> started;

    // Actually add times to the *.series file
    forAll(timeDirs, ti)
    {
        fileName timeDir =
            absolutePath
          ? path.toAbsolute()/timeDirs[ti]
          : relativePath(outputDir, path/timeDirs[ti]);

        wordList files
        (
            readDir
            (
                path/timeDirs[ti],
                fileType::file
            )
        );

        forAll(files, i)
        {
            const fileName& file = files[i];
            const word ext = file.ext();
            if (validExts.found(ext))
            {
                word name(IOobject::member(files[i]));

                Ostream& os = seriesOutputs[file];

                // No comma on the first line
                if (!started.insert(name))
                {
                    os  << ',' << nl;
                }

                // Add time to time series
                os  << indent
                    << '{'
                    << "\"name\":" << timeDir/file << ','
                    << "\"time\":" << times[ti]
                    << '}';
            }
        }
    }

    //- Close brackets in *.series outputs
    forAll(seriesOutputs, i)
    {
        seriesOutputs[i]
            << nl << decrIndent
            << indent << ']' << nl << decrIndent
            << indent << '}' << endl;
    }

    Info<< nl << "Finished" << nl << endl;

    return 0;
}
