/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

#include "vtkTimeSeries.H"
#include "objectRegistry.H"
#include "Time.H"
#include "OFstream.H"
#include "List.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkTimeSeries, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkTimeSeries::vtkTimeSeries(const fileName& path, const label nRemove)
:
    nRemove_(nRemove),
    path_()
{
    fileName casePath(getEnv("FOAM_CASE"));
    wordList caseCmpts(casePath.components());
    wordList cmpts(path.components());
    for (label i = 0; i < cmpts.size() - nRemove_; i++)
    {
        if (caseCmpts[i] != cmpts[i])
        {
            path_ =  path_ / cmpts[i];
        }
    }

    fileNameList dirs(readDir(path_, fileType::directory, true, false));
    forAll(dirs, i)
    {
        IStringStream is(word(dirs[i]));
        token t(is);
        if (t.isNumber())
        {
            this->insert(t.number());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkTimeSeries::~vtkTimeSeries()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::vtkTimeSeries::writeTimeSeries
(
    const fileName& path,
    const fileName& name
) const
{

    OFstream os(path_ / name +".vtk.series");
    wordList cmpts(path.components());
    const_cast<scalarHashSet&>
    (
        dynamicCast<const scalarHashSet&>(*this)
    ).insert
    (
        readScalar((IStringStream(cmpts[cmpts.size() - nRemove_]))())
    );

    // Header
    os
        << '{' << nl
        << "  " << string("file-series-version") << " : " << string("1.0") << ',' << nl
        << "  " << string("files") << " : [";

    scalarList times(this->sortedToc());
    forAll(times, i)
    {
        os  << nl
            << "    { "
            << string("name") << " : "
            << fileName(Foam::name(times[i]) / name + ".vtk") << " ,"
            << string("time") << " : " << times[i]
            << " }";

        if (i != this->size() - 1)
        {
            os << ",";
        }
    }
    os  << nl
        << "  ]" << nl
        << "}" << endl;

    return os.good();
}


// ************************************************************************* //
