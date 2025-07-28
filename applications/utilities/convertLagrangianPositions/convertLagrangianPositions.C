/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020-2021
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
    Utility to convert lagrangian position from baycentric coordinates to
    vector position. The original coordinates are also written so they can
    be reverted if wanted.

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "argList.H"
#include "Time.H"
#include "timeSelector.H"
#include "IFstream.H"
#include "Cloud.H"
#include "particle.H"


namespace Foam
{

enum positionFormat
{
    OLD,
    NEW,
    UNKNOWN
};

positionFormat readCloud
(
    const polyMesh& mesh,
    Cloud<particle>& c,
    IOPosition<Cloud<particle>>& ioP,
    IFstream& is
)
{
    // Start the reading of the list
    // Either read the number or beginning of the list
    token firstToken(is);
    label nParticles = 0;
    if (firstToken.isLabel())
    {
        nParticles = firstToken.labelToken();

        // Read beginning of contents
        is.readBeginList
        (
            "IOPosition<CloudType>::readData(Istream&, CloudType&)"
        );

    }
    else if (firstToken.isPunctuation())
    {
        if (firstToken.pToken() != token::BEGIN_LIST)
        {
            FatalIOErrorInFunction(is)
                << "incorrect first token, '(', found "
                << firstToken.info() << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(is)
            << "incorrect first token, expected <int> or '(', found "
            << firstToken.info() << exit(FatalIOError);
    }

    positionFormat format = UNKNOWN;

    // Read position/coordinates
    if (is.format() == IOstream::ASCII)
    {
        scalarList p0(is);

        // Check the format
        if (p0.size() == 4)
        {
            format = NEW;
        }
        else if (p0.size() == 3)
        {
            format = OLD;
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "Unknown positions format" << endl
                << abort(FatalIOError);
        }

        label pi = 0;
        if (format == NEW)
        {
            token t;
            barycentric p;
            label celli, tetFacei, tetPti;
            while
            (
                is.read(t)
             && !(t.isPunctuation() && t.pToken() == token::END_LIST)
            )
            {
                is.putBack(t);
                if (pi == 0)
                {
                    p[0] = p0[0];
                    p[1] = p0[1];
                    p[2] = p0[2];
                    p[3] = p0[3];
                }
                else
                {
                    is >> p;
                }
                is >> celli >> tetFacei >> tetPti;
                label n = 0;
                c.append
                (
                    new particle(mesh, p, celli, tetFacei, tetPti, n)
                );
                pi++;
            }
        }
        else
        {
            token t;
            vector p;
            label celli;
            while
            (
                is.read(t)
             && !(t.isPunctuation() && t.pToken() == token::END_LIST)
            )
            {
                is.putBack(t);
                if (pi == 0)
                {
                    p[0] = p0[0];
                    p[1] = p0[1];
                    p[2] = p0[2];
                }
                else
                {
                    is >> p;
                }
                is >> celli;
                label n = 0;
                c.append(new particle(mesh, p, celli, n));
                pi++;
            }
        }
    }
    else
    {
        format = NEW;
        for (label i = 0; i < nParticles; i++)
        {
            label n = 0;
            c.append(new particle(is, false));
        }
        // Read beginning of contents
        is.readEndList
        (
            "IOPosition<CloudType>::readData(Istream&, CloudType&)"
        );
    }

    return format;
}

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

    #include "setRegionNames.H"

    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];
        const word& regionDir =
            regionName == polyMesh::defaultRegion
          ? word::null
          : regionName;

        Info<< "\n\nConverting lagrangian positions for region " << regionName << nl
            << endl;

        // Loop over all times
        forAll(timeDirs, timei)
        {
            // Set time for global database
            runTime.setTime(timeDirs[timei], timei);

            Info<< "Time = " << runTime.name() << endl;

            fvMesh mesh
            (
                IOobject
                (
                    regionName,
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                ),
                false
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
            {
                List<fileNameList> procClouds(Pstream::nProcs());
                procClouds[Pstream::myProcNo()] = cloudDirs;
                Pstream::gatherList(procClouds);
                HashSet<fileName> cloudSet;
                forAll(procClouds, proci)
                {
                    cloudSet.insert(procClouds[proci]);
                }
                cloudDirs = cloudSet.toc();
                Pstream::scatter(cloudDirs);
            }

            forAll(cloudDirs, i)
            {
                IOobject positionsIO
                (
                    IOobject
                    (
                        "positions",
                        runTime.name(),
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

                IFstream is(positionsIO.objectPath(false));
                bool write = true;
                positionFormat format = NEW;
                if (is.good())
                {
                    positionsIO.readHeader(is);
                    cloudType = positionsIO.headerClassName();

                    format = readCloud(mesh, c, ioP, is);
                }
                else
                {
                    write = false;
                }

                write =
                    write
                 && (
                        (!revert && (format == NEW))
                     || (revert && (format == OLD))
                    );

                if (write)
                {
                    Info << "\tWriting positions file" << endl;

                    OFstream positionsOS(positionsIO.objectPath(false));
                    positionsIO.writeHeader
                    (
                        positionsOS,
                        revert
                      ? cloudType.replaceAll("Cloud", '\0')
                      : cloudType + "Cloud"
                    );

                    positionsOS  << c.size() << nl << token::BEGIN_LIST << nl;
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
                                << pIter().position(mesh)
                                << token::SPACE << pIter().cell()
                                << nl;
                        }
                    }
                    positionsOS  << token::END_LIST << endl;
                }
                else
                {
                    Info<< "\tNot writing positions file. Already converted" << endl;
                }
                Info<< endl;
            }
        }
    }
}
