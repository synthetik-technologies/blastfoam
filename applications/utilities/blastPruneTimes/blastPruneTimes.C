/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
-------------------------------------------------------------------------------
03-02-2023:                 | Modified foamListTimes to allow for selecting
    Synthetik Applied       | times within an interval and with prescribed
    Technologie             | spacing
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    blastPruneTimes

Description
    Select time between a start and end time with a given spacing

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;
int main(int argc, char *argv[])
{
    writeInfoHeader = false;

    argList::addNote("List times from start : spacing : end");
    timeSelector::addOptions(true, true);
    argList::noParallel();
    argList::addBoolOption
    (
        "processor",
        "list times from processor0/ directory"
    );
    argList::addBoolOption
    (
        "rm",
        "remove selected time directories"
    );
    argList::addBoolOption
    (
        "withFunctionObjects",
        "execute functionObjects"
    );

    argList::addOption
    (
        "startTime",
        "Start of interval"
    );
    argList::addOption
    (
        "endTime",
        "End of interval"
    );
    argList::addOption
    (
        "dt",
        "Time spacing"
    );
    argList::addOption
    (
        "n",
        "Number of times to save"
    );

    argList::addOption
    (
        "tolerance",
        "Tolerance for finding times"
    );

    argList::addBoolOption
    (
        "invert",
        "Select non-interval times"
    );

    argList::addBoolOption
    (
        "Clean",
        "Remove times not in the processor0 folder"
    );

    argList::addBoolOption("v", "Verbose");


    #include "addRegionOption.H"
    #include "setRootCase.H"

    const bool verbose = args.optionFound("v");

    label nProcs = 1;
    PtrList<Time> databases(1);
    if (args.optionFound("processor"))
    {
        nProcs = fileHandler().nProcs(args.path());
        if (!nProcs)
        {
            FatalErrorInFunction
                << "Trying to show time for a parallel case, but the case" << nl
                << "has not been deomcomposed." << endl
                << abort(FatalError);
        }
        databases.setSize(nProcs);
        forAll(databases, proci)
        {
            databases.set
            (
                proci,
                new Time
                (
                    Time::controlDictName,
                    args.rootPath(),
                    args.caseName()/fileName("processor" + Foam::name(proci))
                )
            );
        }
    }
    else
    {
        databases.set
        (
            0,
            new Time
            (
                Time::controlDictName,
                args.rootPath(),
                args.caseName()
            )
        );
    }

    instantList times = timeSelector::select
    (
        databases[0].times(),
        args
    );

    scalar minDt = great;
    for (label i = 1; i < times.size(); i++)
    {
        minDt = min(minDt, times[i].value() - times[i-1].value());
    }

    DynamicList<instant> timesToKeep(times.size());
    DynamicList<instant> prunedTimes(times.size());
    const scalar startTime =
        args.optionLookupOrDefault("startTime", databases[0].beginTime().value());
    const scalar endTime =
        args.optionLookupOrDefault("endTime", databases[0].endTime().value());
    scalar dt = endTime - startTime;
    if (args.optionFound("dt"))
    {
        dt = args.optionRead<scalar>("dt");
    }
    else if (args.optionFound("n"))
    {
        dt /= args.optionRead<scalar>("n");
    }

    const scalar tolerance = args.optionLookupOrDefault
    (
        "tolerance",
        minDt*1e-3
    );

    scalar nextTime = startTime;

//     // Make sure the first time in the selected times is less than the starting time
//     if (times.size())
//     {
//         scalar tByDt(times[0].value()/dt);
//         if (mag(label(tByDt) - tByDt) < tolerance)
//         {
//             nextTime = times[0].value();
//         }
//         else
//         {
//             forAll(times, ti)
//             {
//                 if (times[ti].value() > startTime)
//                 {
//                     nextTime += dt;
//                 }
//                 else
//                 {
//                     break;
//                 }
//             }
//         }
//     }

    // Space times
    forAll(times, ti)
    {
        if (times[ti].value() > endTime+tolerance)
        {
            break;
        }
        else if (times[ti].value() < startTime-tolerance)
        {
            continue;
        }

        scalar tByDt(times[ti].value()/dt);
        if (mag(round(tByDt) - tByDt) < tolerance)
        {
            timesToKeep.append(times[ti]);
        }
        else
        {
            prunedTimes.append(times[ti]);
        }
    }

    const List<instant>& selectedTimes =
        args.optionFound("invert")
      ? prunedTimes
      : timesToKeep;

    if (args.optionFound("rm"))
    {
        if (args.optionFound("processor"))
        {
            HashSet<word> masterTimes;
            forAll(times, ti)
            {
                masterTimes.insert(times[ti].name());
            }

            for (label proci=0; proci<nProcs; proci++)
            {
                const fileName procPath
                (
                    args.path()/(word("processor") + name(proci))
                );

                forAll(selectedTimes, ti)
                {
                    const fileName procTimePath
                    (
                        fileHandler().filePath(procPath/selectedTimes[ti].name())
                    );
                    if (isDir(procTimePath))
                    {
                        if (verbose)
                        {
                            Info<< "Removing " << procTimePath << endl;
                        }
                        rmDir(procTimePath);
                    }
                }

                instantList procTimes
                (
                    timeSelector::select
                    (
                        databases[proci].times(),
                        args
                    )
                );
                forAll(procTimes, ti)
                {
                    const fileName procTimePath
                    (
                        fileHandler().filePath
                        (
                            procPath/procTimes[ti].name()
                        )
                    );
                    if
                    (
                        !masterTimes.found(procTimes[ti].name())
                     && isDir(procTimePath)
                    )
                    {
                        if (verbose)
                        {
                            Info<< "removing " << procTimePath << endl;
                        }
                        rmDir(procTimePath);
                    }
                }
            }
        }
        else
        {
            forAll(selectedTimes, ti)
            {
                const fileName timePath
                (
                    fileHandler().filePath(args.path()/selectedTimes[ti].name())
                );

                if (verbose)
                {
                    Info<< "Removing " << timePath << endl;
                }
                rmDir(timePath);
            }
        }
    }
    else
    {
        forAll(selectedTimes, ti)
        {
            Info<< selectedTimes[ti].name() << endl;
        }
    }

    return 0;

}


// ************************************************************************* //
