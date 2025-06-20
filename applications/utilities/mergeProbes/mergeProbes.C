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
    Utility to merge probe files from multiple start times

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Istream.H"
#include "IFstream.H"
#include "OFstream.H"
#include "SortableList.H"
#include "scalarList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Merges probe file started from different times\n\n"
    );

    argList::validArgs.append("probeDir");
    argList::addBoolOption
    (
        "force",
        "Remove VTK directory if currently present"
    );
    argList::addOption
    (
        "fieldNames",
        "Name of probe to merge"
    );
    argList::addOption
    (
        "times",
        "List of times to include"
    );
    argList::addBoolOption
    (
        "local",
        "Path to probe is directly specified"
    );

    #include "setRootCase.H"

    bool force(args.optionFound("force"));
    wordRe probeDirName(args.argRead<wordRe>(1), wordRe::compOption::detect);
    wordList fieldNames(args.optionLookupOrDefault("fieldNames", wordList()));

    // Create the probe databases
    fileName postProcessingDir
    (
        args.optionFound("local")
      ? probeDirName
      : args.rootPath()/args.caseName()/fileName("postProcessing")
    );

    //- Read a directory and return the entries as a string list
    fileNameList ppDirs
    (
        readDir(postProcessingDir, fileType::directory, false, false)
    );
    fileNameList probeDirs;
    forAll(ppDirs, i)
    {
        if (probeDirName.match(ppDirs[i].name()))
        {
            probeDirs.append(postProcessingDir/ppDirs[i]);
        }
    }

    forAll(probeDirs, i)
    {
        const fileName& probeDir = probeDirs[i];
        if (!isDir(probeDir))
        {
            FatalErrorInFunction
                << "Provided probe directory, " << probeDir
                << ", does not exist" << endl
                << exit(FatalError);
        }


        Info<< "Merge probes in directory: " << probeDir.name()
            << incrIndent << endl;

        wordList times;
        SortableList<scalar> sTimes;
        if (args.optionFound("times"))
        {
            sTimes = args.optionRead<scalarList>("times");
            sTimes.sort();

            times.setSize(sTimes.size());
            forAll(sTimes, ti)
            {
                times[ti] = Foam::Time::timeName(sTimes[ti]);
            }
        }
        else
        {
            times = wordList(readDir(probeDir, fileType::directory));
            if (!times.size())
            {
                FatalErrorInFunction
                    << "No times were found for probe " << probeDir << endl
                    << exit(FatalError);
            }
            sTimes.setSize(times.size());

            // Sort times
            forAll(sTimes, ti)
            {
                IStringStream is(times[ti]);
                sTimes[ti] = readScalar(is);
            }
            sTimes.sort();
            wordList oldTimes(times);
            forAll(sTimes, ti)
            {
                times[ti] = oldTimes[sTimes.indices()[ti]];
            }
        }
        Info<< indent << "Merging times:" << nl << incrIndent;
        forAll(times, ti)
        {
            Info << indent << times[ti] << nl;
        }
        Info<< decrIndent << endl;
        sTimes.append(great);

        // Get full list of probes
        if (!args.optionFound("fieldNames"))
        {
            fileName probe0Dir(probeDir/times[0]);
            fieldNames = wordList(readDir(probe0Dir, fileType::file));
        }
        if (!force)
        {
            wordList writtenProbes;
            forAll(fieldNames, probei)
            {
                if (!isFile(probeDir/fieldNames[probei]))
                {
                    writtenProbes.append(fieldNames[probei]);
                }
                else
                {
                    Warning << nl
                        << (probeDir/fieldNames[probei])
                        << " already found. Skipping probe." << nl
                        << endl;
                }
            }
            fieldNames = writtenProbes;
        }

        Info<< indent << "Probes to merge:" << nl << incrIndent;
        forAll(fieldNames, probei)
        {
            Info << indent << fieldNames[probei] << nl;
        }
        Info<< decrIndent << endl;


        // Create outputs
        PtrList<OFstream> outputs(fieldNames.size());
        forAll(outputs, probei)
        {
            outputs.set(probei, new OFstream(probeDir/fieldNames[probei]));
        }

        scalar nextTime = -1.0;
        bool header = true;
        forAll(times, timei)
        {
            nextTime = sTimes[timei + 1];
            fileName probeTimeDir(probeDir/times[timei]);

            forAll(fieldNames, probei)
            {
                IFstream stream(probeTimeDir/fieldNames[probei]);

                while (stream.good())
                {
                    string line;
                    stream.getLine(line);

                    if (line[0] == '#')
                    {
                        if (header)
                        {
                            outputs[probei] << word(line) << nl;
                        }
                        continue;
                    }
                    header = false;

                    IStringStream is(line);
                    scalar t = readScalar(is);

                    if (t < nextTime)
                    {
                        outputs[probei] << word(line) << nl;
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        Info<< decrIndent << endl;
    }

    Info<< nl << "Done." << endl;
}
