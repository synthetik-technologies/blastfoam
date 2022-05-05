/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 Synthetik Applied Technologies
     \\/     M anipulation  |
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

#include "PtrListDictionary.H"
#include "OFstream.H"
#include "argList.H"
#include "Time.H"
#include "timeSelector.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    timeSelector::addOptions(true, false);

    argList::addBoolOption
    (
        "relative",
        "Create the VTK series using relative file paths."
    );

    #include "setRootCase.H"
    #include "createTime.H"

    //- Select time
    runTime.functionObjects().off();
    instantList timeDirs = timeSelector::selectIfPresent(runTime, args);

    bool relativePath = args.optionFound("relative");

    // Create the processor databases
    fileName postProcessDir
    (
        args.caseName()/fileName(word("postProcessing"))
    );
    wordList dirs
    (
        readDir
        (
            args.rootPath()/fileName(postProcessDir),
            fileType::directory
        )
    );

    fileName VTKDir(args.rootPath()/args.caseName()/fileName("VTK"));
    mkDir(VTKDir);

    forAll(dirs, diri)
    {
        Time database
        (
            Time::controlDictName,
            args.rootPath(),
            postProcessDir/fileName(dirs[diri]),
            word("../../system")
        );

        // Use the times list from the master processor
        // and select a subset based on the command-line options
        instantList timeDirs = timeSelector::select
        (
            database.times(),
            args
        );

        //- List of vtk.series files
        PtrListDictionary<OFstream> seriesOutputs(0);

        //- Check all directories for vtk files and create series outputs
        forAll(timeDirs, timeI)
        {
            fileName timeDir
            (
                args.rootPath()
               /postProcessDir
               /fileName(dirs[diri])
               /fileName(timeDirs[timeI].name())
            );
            wordList files
            (
                readDir
                (
                    timeDir,
                    fileType::file
                )
            );

            forAll(files, i)
            {
                if (IOobject::group(files[i]) == "vtk")
                {
                    word name(IOobject::member(files[i]));

                    if (!seriesOutputs.found(name))
                    {
                        Info<< "Making time series file, " << name << ".vtk.series"
                            << " in " << word(VTKDir) << endl;

                        seriesOutputs.resize(seriesOutputs.size() + 1);
                        seriesOutputs.set
                        (
                            seriesOutputs.size() - 1,
                            name,
                            new OFstream
                            (
                                VTKDir
                               /name + ".vtk.series"
                            )
                        );

                        // Header
                        seriesOutputs[name]
                            << '{' << nl
                            << "  \"file-series-version\" : \"1.0\"," << nl
                            << "  \"files\" : [";
                    }
                }
            }
        }


        HashSet<word> start;

        forAll(timeDirs, timeI)
        {
            fileName timeDir
            (
                args.rootPath()
               /postProcessDir
               /fileName(dirs[diri])
               /fileName(timeDirs[timeI].name())
            );
            fileName relDir
            (
                fileName("../postProcessing")
               /fileName(dirs[diri])
               /fileName(timeDirs[timeI].name())
            );
            wordList files
            (
                readDir
                (
                    timeDir,
                    fileType::file
                )
            );

            forAll(files, i)
            {
                if (IOobject::group(files[i]) == "vtk")
                {
                    word name(IOobject::member(files[i]));

                    // No comma on the first line
                    if (start.found(name))
                    {
                        seriesOutputs[name]
                            << ",";
                    }
                    else
                    {
                        start.set(name);
                    }
                    // Add time to time series
                    seriesOutputs[name]
                        << nl
                        << "    { "
                        << "\"name\" : "
                        << "\"" <<
                            (relativePath ?
                            word(relDir/fileName(files[i])) :
                            word(timeDir/fileName(files[i])))
                        << "\" ,"
                        << "\"time\" : "
                        << timeDirs[timeI].value()
                        << " }";
                }
            }
        }

        //- Close brackets in all series outputs
        forAll(seriesOutputs, i)
        {
            seriesOutputs[i]
                << nl
                << "  ]" << nl
                << "}" << endl;
        }
    }
}
