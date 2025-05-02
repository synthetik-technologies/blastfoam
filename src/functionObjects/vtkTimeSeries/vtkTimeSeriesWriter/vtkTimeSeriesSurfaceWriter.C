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
#include "vtkWritePolyData.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(vtkTimeSeriesSurfaceWriter, 0);
    addToRunTimeSelectionTable(surfaceWriter, vtkTimeSeriesSurfaceWriter, word);
    addToRunTimeSelectionTable(surfaceWriter, vtkTimeSeriesSurfaceWriter, dict);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkTimeSeriesSurfaceWriter::~vtkTimeSeriesSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtkTimeSeriesSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const wordList& fieldNames,
    const bool writePointValues
    #define FieldTypeValuesConstArg(Type, nullArg) \
        , const UPtrList<const Field<Type>>& field##Type##Values
    FOR_ALL_FIELD_TYPES(FieldTypeValuesConstArg)
    #undef FieldTypeValuesConstArg
) const
{
    const fileName surfaceDir(outputDir/surfaceName);

    if (!isDir(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    if (!timeSeries_.valid())
    {
        timeSeries_.set(new vtkTimeSeries(outputDir, 1, true)); // Read
    }

    vtkWritePolyData::write
    (
        outputDir/surfaceName + ".vtk",
        "sampleSurface",
        writeFormat_ == IOstream::BINARY,
        points,
        labelList(),
        edgeList(),
        faces,
        fieldNames,
        boolList(fieldNames.size(), writePointValues),
        UPtrList<const Field<label>>(fieldNames.size())
        #define FieldTypeValuesParameter(Type, nullArg) , field##Type##Values
        FOR_ALL_FIELD_TYPES(FieldTypeValuesParameter)
        #undef FieldTypeValuesParameter
    );

    // Remove 1 level from the end of the path
    timeSeries_->insertFromPath(outputDir, 1);
    timeSeries_->writeTimeSeries(surfaceName);
}

// ************************************************************************* //
