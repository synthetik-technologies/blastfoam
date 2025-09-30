/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "blastSampledSurfaces.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ListListOps.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::functionObjects::blastSampledSurfaces::writeSurface
(
    const Field<Type>& values,
    const label surfI,
    const word& fieldName,
    const fileName& outputDir,
    const bool isNodeValue
)
{
    const blastSampledSurface& s = operator[](surfI);

    if (Pstream::parRun())
    {
        // Collect values from all processors
        List<Field<Type>> gatheredValues(Pstream::nProcs());
        gatheredValues[Pstream::myProcNo()] = values;
        Pstream::gatherList(gatheredValues);

        if (Pstream::master())
        {
            // Combine values into single field
            Field<Type> allValues
            (
                ListListOps::combine<Field<Type>>
                (
                    gatheredValues,
                    accessOp<Field<Type>>()
                )
            );

            // Renumber (point data) to correspond to merged points
            if (mergeList_[surfI].pointsMap.size() == allValues.size())
            {
                inplaceReorder(mergeList_[surfI].pointsMap, allValues);
                allValues.setSize(mergeList_[surfI].points.size());
            }

            // Write to time directory under outputPath_
            // skip surface without faces (eg, a failed cut-plane)
            if (mergeList_[surfI].faces.size())
            {
                formatter_->write
                (
                    outputDir,
                    s.name(),
                    mergeList_[surfI].points,
                    mergeList_[surfI].faces,
                    fieldName,
                    allValues,
                    isNodeValue
                );
            }
        }
    }
    else
    {
        // Write to time directory under outputPath_
        // skip surface without faces (eg, a failed cut-plane)
        if (s.faces().size())
        {
            formatter_->write
            (
                outputDir,
                s.name(),
                s.points(),
                s.faces(),
                fieldName,
                values,
                isNodeValue
            );
        }
    }
}


template<class Type>
void Foam::functionObjects::blastSampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    // interpolator for this field
    autoPtr<interpolation<Type>> interpolatorPtr;

    const word& fieldName = vField.name();
    const fileName outputDir = outputPath_/vField.time().timeName();

    forAll(*this, surfI)
    {
        const blastSampledSurface& s = operator[](surfI);

        Field<Type> values;

        if (s.interpolate())
        {
            if (interpolatorPtr.empty())
            {
                interpolatorPtr = interpolation<Type>::New
                (
                    interpolationScheme_,
                    vField
                );
            }

            values = s.interpolate(interpolatorPtr());
        }
        else
        {
            values = s.sample(vField);
        }

        writeSurface<Type>(values, surfI, fieldName, outputDir, s.interpolate());
    }
}


template<class Type>
void Foam::functionObjects::blastSampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& sField
)
{
    const word& fieldName   = sField.name();
    const fileName outputDir = outputPath_/sField.time().timeName();

    forAll(*this, surfI)
    {
        const blastSampledSurface& s = operator[](surfI);
        Field<Type> values(s.sample(sField));
        writeSurface<Type>(values, surfI, fieldName, outputDir, s.interpolate());
    }
}

template<class Type>
void Foam::functionObjects::blastSampledSurfaces::sampleAndWrite
(
    const GeometricField<Type, pointPatchField, pointMesh>& pField
)
{
    const word& fieldName = pField.name();
    const fileName outputDir = outputPath_/pField.time().timeName();

    forAll(*this, surfI)
    {
        const blastSampledSurface& s = operator[](surfI);
        Field<Type> values(s.sample(pField));
        writeSurface<Type>(values, surfI, fieldName, outputDir, true);
    }
}


template<class GeoField>
void Foam::functionObjects::blastSampledSurfaces::sampleAndWrite
(
    const IOobjectList& objects
)
{
    HashTable<const GeoField*> fields(mesh_.lookupClass<GeoField>());
    forAllConstIter(typename HashTable<const GeoField*>, fields, iter)
    {
        const word& fieldName = iter.key();

        if (findIndex(fieldSelection_, fieldName) >= 0)
        {
            if (verbose_)
            {
                Info<< "sampleAndWrite: " << iter()->name() << endl;
            }
            sampleAndWrite(*(iter()));
        }
    }
}

// ************************************************************************* //
