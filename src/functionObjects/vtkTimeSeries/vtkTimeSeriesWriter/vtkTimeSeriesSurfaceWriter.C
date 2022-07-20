/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2022-05-09: Synthetik Applpied Technologies : Added support for vtkTimeSeries
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "vtkTimeSeriesSurfaceWriter.H"
#include "OFstream.H"
#include "boolList.H"
#include "OSspecific.H"
#include "makeSurfaceWriterMethods.H"
#include "vtkWritePolyData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(vtkTimeSeriesSurfaceWriter);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::vtkTimeSeriesSurfaceWriter::Write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    if (!timeSeries_.valid())
    {
        timeSeries_.set(new vtkTimeSeries(outputDir, 1));
    }

    vtkWritePolyData::write
    (
        outputDir/fieldName + '_' + surfaceName + ".vtk",
        "sampleSurface",
        writeFormat_ == IOstream::BINARY,
        points,
        labelList(),
        edgeList(),
        faces,
        fieldName,
        isNodeValues,
        values
    );
    timeSeries_->writeTimeSeries(outputDir, fieldName + '_' + surfaceName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkTimeSeriesSurfaceWriter::vtkTimeSeriesSurfaceWriter
(
    const IOstream::streamFormat writeFormat
)
:
    surfaceWriter(writeFormat)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkTimeSeriesSurfaceWriter::~vtkTimeSeriesSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkTimeSeriesSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    vtkWritePolyData::write
    (
        outputDir/surfaceName + ".vtk",
        "sampleSurface",
        writeFormat_ == IOstream::BINARY,
        points,
        labelList(),
        edgeList(),
        faces
    );
}


// Create write methods
defineSurfaceWriterWriteFields(Foam::vtkTimeSeriesSurfaceWriter);


// ************************************************************************* //
