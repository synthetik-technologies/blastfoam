/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2022
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
    Utility to calculate the impulse given a pressure probe

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Istream.H"
#include "IFstream.H"
#include "OFstream.H"
#include "SortableList.H"

using namespace Foam;

fileName findProbeDir(const fileName& dir, const word& pName)
{
    if (isFile(dir/pName))
    {
        return dir;
    }


    fileNameList timeDirs(readDir(dir, fileType::directory));

    scalar minTime = great;
    fileName minTimeName;
    forAll(timeDirs, ti)
    {
        token t((IStringStream(timeDirs[ti]))());
        if (t.isNumber() && isFile(dir/timeDirs[ti]/pName))
        {
            if (t.number() < minTime)
            {
                minTime = t.number();
                minTimeName = timeDirs[ti];
            }
        }
    }

    if (!minTimeName.empty())
    {
        return dir/minTimeName;
    }
    return fileName::null;
}

fileNameList collectFiles(const fileName& dir)
{
    HashSet<fileName> collectedFiles
    (
        readDir(dir, fileType::file)
    );

    fileNameList timeDirs(readDir(dir, fileType::directory));

    forAll(timeDirs, ti)
    {
        token t((IStringStream(timeDirs[ti]))());
        if (t.isNumber())
        {
            collectedFiles.insert
            (
                readDir(dir/timeDirs[ti], fileType::file)
            );
        }
    }

    return collectedFiles.toc();
}


int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote
    (
        "Calculates impulse based on pressure from probe files\n\n"
    );

    argList::validArgs.append("probeDir");
    argList::addOption
    (
        "p",
        "Name of pressure file"
    );
    argList::addOption
    (
        "impulse",
        "Name of impulse file"
    );
    argList::addOption
    (
        "pRef",
        "Reference pressure [Pa]"
    );
    argList::addBoolOption
    (
        "force",
        "overwrite exisiting impulse file"
    );

    #include "setRootCase.H"

    fileName probeDirName(args.argRead<fileName>(1));
    fileName probeDir("postProcessing"/probeDirName);

    // Default pressure name
    fileName pName;
    if (args.optionFound("p"))
    {
        pName = args.optionRead<fileName>("p");
        fileName newProbeDir = findProbeDir(probeDir, pName);

        if (newProbeDir.empty())
        {
            // Read all available probed fields
            fileNameList probeFields(collectFiles(probeDir));

            FatalErrorInFunction
                << pName << " was not found in " << probeDir << nl
                << "valid probe fields are:" << nl
                << probeFields << endl
                << abort(FatalError);
        }
        probeDir = newProbeDir;
    }
    else
    {
        fileName newProbeDir = findProbeDir(probeDir, "overpressure");
        if (!newProbeDir.empty())
        {
            pName = "overpressure";
            probeDir = newProbeDir;
        }
        else
        {
            newProbeDir = findProbeDir(probeDir, "p");
            if (!newProbeDir.empty())
            {
                pName = "p";
                probeDir = newProbeDir;
            }
            else
            {
                FatalErrorInFunction
                    << "Could not detemine a pressure field to calulate impulse with." << nl
                    << "please specify a field using \'-p pName\'" << endl
                    << abort(FatalError);
            }
        }
    }
    Info<< "Using " << string(pName) << " in " << probeDir << endl;

    fileName impulseName
    (
        args.optionLookupOrDefault
        (
            "impulse",
            IOobject::groupName("impulse", IOobject::group(pName))
        )
    );

    // Check if the file exists
    if (isFile(probeDir/impulseName))
    {
        // Overwrite
        if (args.optionFound("force"))
        {
            Info<< "Overwriting " << (probeDir/impulseName) << endl;
        }

        // Return since not overwriting
        else
        {
            WarningInFunction
                << (probeDir/impulseName) << " already exisits." << nl
                << "    Use \'-force\' to overwrite" << nl << endl;

            return 0;
        }
    }
    else
    {
        Info<< "Writing " << impulseName << " to " << probeDir << endl;
    }

    // Open ofstream
    OFstream impulseStream(probeDir/impulseName);

    // Read the pressure probes and check if it exists
    IFstream pstream(probeDir/pName);
    if (!pstream.good())
    {
        FatalIOErrorInFunction(pstream)
            << "Could not find " << probeDir/pName << endl
            << abort(FatalIOError);
    }


    string line;

    // Read the header and determine the number of probes
    label nProbes = 0;
    pstream.getLine(line);
    while (line[0] == '#')
    {
        nProbes++;
        impulseStream << word(line) << nl;
        pstream.getLine(line);
    }

    // Last commented line is
    // # time p0 p1 p2 ...
    nProbes--;

    // Read the first line
    IStringStream is(line);

    // Read the time
    scalar t(readScalar(is));

    // Read the pressures
    scalarField p(nProbes);
    forAll(p, i)
    {
        is >> p[i];
    }

    // Old values
    scalar tOld(t);
    scalarField pOld(p);
    scalarField impulse(nProbes, 0.0);

    // Set reference pressure
    scalar pRef;
    if (args.optionFound("pRef"))
    {
        pRef = args.optionRead<scalar>("pRef");
    }
    else
    {
        scalar pAvg = 0;
        forAll(pOld, i)
        {
            pAvg += pOld[i];
        }
        pRef = pAvg/scalar(pOld.size());
    }
    Info<< "Using reference pressure of " << pRef << " Pa" << endl;

    // Initial time
    scalar t0 = t;

    while (pstream.good())
    {
        pstream.getLine(line);
        IStringStream is(line);

        // Save old time and pressure
        tOld = t;
        pOld = p;

        // Read next time
        is >> t;

        // Skip dts less than small
        if (mag(t - tOld) < small)
        {
            continue;
        }
        scalar dt = t - tOld;

        // Read pressures
        forAll(p, i)
        {
            is >> p[i];
        }

        // Update impulse with average pressure
        impulse += (0.5*(p + pOld) - pRef)*dt;

        // Write to the output
        impulseStream << t << " ";
        forAll(p, i)
        {
            impulseStream<< impulse[i] << " ";
        }
        impulseStream << endl;
    }

    Info<< nl << "Integrated from "
        << t0 << " s to " << t << " s" << nl
        << nl << "Done." << endl;
}
