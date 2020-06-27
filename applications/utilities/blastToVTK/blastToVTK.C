/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
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
    Utility to create symbolic link to post processing vtk files for easier
    viewing in paraview

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "force",
        "Remove VTK directory if currently present"
    );
    argList::addBoolOption
    (
        "useTimeName",
        "Rename VTK files with time step name"
    );
    argList::addBoolOption
    (
        "hardCopy",
        "Hard copy of VTK files"
    );

    #include "setRootCase.H"

    bool useTime(args.optionFound("useTimeName"));
    bool hardCopy(args.optionFound("hardCopy"));

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
    if (args.optionFound("force"))
    {
        Info<< "Remving old VTK directory" << endl;
        rmDir(VTKDir);
    }
    mkDir(VTKDir);

    forAll(dirs, diri)
    {
        Info<< "Linking vtk files from " << dirs[diri] << " directory"
            << " to " << word(VTKDir) << endl;

        mkDir(VTKDir/fileName(dirs[diri]));

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
                    word index;
                    if (useTime)
                    {
                        index = timeDirs[timeI].name();
                    }
                    else
                    {
                        index = Foam::name(timeI);
                        while (index.size() < 6)
                        {
                            index = '0' + index;
                        }
                    }

                    word name(IOobject::member(files[i]));
                    if (hardCopy)
                    {
                        cp
                        (
                            timeDir/fileName(files[i]),
                            VTKDir
                           /fileName(dirs[diri])
                           /fileName
                            (
                                name
                              + word("_")
                              + index
                              + word(".vtk")
                            )
                        );
                    }
                    else
                    {
                        ln
                        (
                            timeDir/fileName(files[i]),
                            VTKDir
                           /fileName(dirs[diri])
                           /fileName
                            (
                                name
                              + word("_")
                              + index
                              + word(".vtk")
                            )
                        );
                    }
                }
            }
        }
    }
}
