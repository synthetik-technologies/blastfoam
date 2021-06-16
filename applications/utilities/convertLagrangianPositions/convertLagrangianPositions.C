/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
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
    Utility to convert lagrangian position from baycentric coordinates to
    vector position. The original coordinates are also written so they can
    be reverted if wanted.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "argList.H"
#include "Time.H"
#include "regionProperties.H"
#include "timeSelector.H"
#include "Cloud.H"
#include "particle.H"


namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<particle>, 0);
}



using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert lagrangian position file to be read in Paraview"
    );

    // Enable -constant ... if someone really wants it
    // Enable -withZero to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);

    argList::addBoolOption("revert", "Revert from positions to coordinates");

    #include "addRegionOption.H"
    #include "addAllRegionsOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    // Use the times list from the master processor
    // and select a subset based on the command-line options
    instantList timeDirs = timeSelector::select
    (
        runTime.times(),
        args
    );

    if (timeDirs.empty())
    {
        WarningInFunction << "No times selected" << endl;
        exit(1);
    }

    bool revert = args.optionFound("revert");

    const wordList regionNames(selectRegionNames(args, runTime));

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir = Foam::regionDir(regionName);

        Info<< "\n\nConverting lagrangian positions for region " << regionName << nl
            << endl;

        // Loop over all times
        forAll(timeDirs, timei)
        {
            // Set time for global database
            runTime.setTime(timeDirs[timei], timei);

            Info<< "Time = " << runTime.timeName() << endl << endl;

            fvMesh mesh
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            );
            mesh.tetBasePtIs();

            fileName lagrangianDir
            (
                fileHandler().filePath
                (
                    runTime.timePath()
                  / regionDir
                  / cloud::prefix
                )
            );

            fileNameList cloudDirs;
            if (!lagrangianDir.empty())
            {
                cloudDirs = fileHandler().readDir
                (
                    lagrangianDir,
                    fileType::directory
                );
            }

            forAll(cloudDirs, i)
            {

                IOobject positionsIO
                (
                    IOobject
                    (
                        "positions",
                        runTime.timeName(),
                        cloud::prefix/cloudDirs[i],
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                );
                IOobject coordinatesIO
                (
                    IOobject
                    (
                        "coordinates",
                        runTime.timeName(),
                        cloud::prefix/cloudDirs[i],
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                );
                word cloudType;
                IDLList<particle> tmp;
                Cloud<particle> c(mesh, cloudDirs[i], tmp);
                IOPosition<Cloud<particle>> ioP(c);
                {
                    IFstream pis(positionsIO.objectPath());
                    positionsIO.readHeader(pis);

                    IFstream cis(coordinatesIO.objectPath());
                    coordinatesIO.readHeader(cis);

                    cloudType = positionsIO.headerClassName();
                    ioP.headerClassName() = cloudType;

                    if (revert)
                    {
                        ioP.readData(cis, c);
                    }
                    else
                    {
                        ioP.readData(pis, c);
                    }
                }

                OFstream positionsOS(positionsIO.objectPath());
                positionsIO.writeHeader(positionsOS, cloudType);
                OFstream coordinatesOS(coordinatesIO.objectPath());
                coordinatesIO.writeHeader(coordinatesOS, cloudType);

                positionsOS  << c.size() << nl << token::BEGIN_LIST << nl;
                coordinatesOS  << c.size() << nl << token::BEGIN_LIST << nl;
                forAllConstIter(Cloud<particle>, c, pIter)
                {
                    if (revert)
                    {
                        positionsOS
                            << pIter().coordinates()
                            << token::SPACE << pIter().cell()
                            << token::SPACE << pIter().tetFace()
                            << token::SPACE << pIter().tetPt()
                            << nl;
                    }
                    else
                    {
                        positionsOS
                            << pIter().position()
                            << token::SPACE << pIter().cell()
                            << nl;

                        coordinatesOS
                            << pIter().coordinates()
                            << token::SPACE << pIter().cell()
                            << token::SPACE << pIter().tetFace()
                            << token::SPACE << pIter().tetPt()
                            << nl;
                    }
                }
                positionsOS  << token::END_LIST << endl;
                coordinatesOS  << token::END_LIST << endl;

            }
        }
    }
}
