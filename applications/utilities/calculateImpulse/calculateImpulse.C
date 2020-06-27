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
    Utility to calculate the impulse given a pressure probe

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Istream.H"
#include "IFstream.H"
#include "OFstream.H"
#include "SortableList.H"

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "probeDir",
        "Name of probe directory"
    );
    argList::addOption
    (
        "name",
        "Name of pressure field"
    );
    argList::addOption
    (
        "pRef",
        "Reference pressure [Pa]"
    );

    #include "setRootCase.H"

    word name(args.optionLookupOrDefault("name", word("p")));
    word probeName(args.option("probeDir"));
    IStringStream pRefStream(IStringStream(args.option("pRef")));
    scalar pRef(readScalar(pRefStream));
    fileName probeDir
    (
        args.rootPath()/args.caseName()/fileName(word("postProcessing"))/probeName
    );

    fileName pFile(probeDir);
    if (!isFile(probeDir/name))
    {
        if (isDir(probeDir/"0"))
        {
            probeDir = probeDir/"0";
        }
    }

    pFile = probeDir/name;
    if (!isFile(pFile))
    {
        FatalErrorInFunction
            << name << " Was not found in " << probeDir << "."
            << abort(FatalError);
    }
    Info<<pFile<<endl;

    bool header = true;
    scalar t = 0;
    scalar tOld = 0;
    scalarField p;
    scalarField pOld;

    scalarField impulse;
    OFstream impulseStream(probeDir/"impulse");

    // Get number of probes
    label nProbes = 0;
    {
        IFstream pstream(pFile);
        while (pstream.good())
        {
            string line;
            pstream.getLine(line);

            if (line[0] == '#')
            {
                if (header)
                {
                    impulseStream << word(line) << nl;
                }
                continue;
            }
            if (header)
            {
                IStringStream is(line);
                while (is.good())
                {
                    readScalar(is);
                    nProbes++;
                }
                pstream.rewind();
                nProbes--;
                header = false;
                p.resize(nProbes);
                p = pRef;
                pOld.resize(nProbes);
                impulse.resize(nProbes);
                impulse = 0;
            }

            IStringStream is(line);
            tOld = t;
            is >> t;
            scalar dt = t - tOld;

            pOld = p;
            label i = 0;
            while (is.good())
            {
                is >> p[i];
                i++;
            }
            impulse += (0.5*(p + pOld) - pRef)*dt;
            impulseStream << t << " ";
            for (label i = 0; i  < nProbes; i++)
            {
                impulseStream<< impulse[i] << " ";
            }
            impulseStream << endl;
        }
    }

    Info<< nl << "Done." << endl;
}
