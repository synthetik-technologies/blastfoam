/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
13-05-2020 Jeff Heylmun:    | Store field maxes over all times
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

#include "fieldMax.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMax, 0);
    addToRunTimeSelectionTable(functionObject, fieldMax, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMax::fieldMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    fieldNames_(dict.lookup("fields")),
    maxFieldNames_(fieldNames_.size()),

    cellMap_(nullptr),
    rCellMap_(nullptr),

    vScalarFields_(0),
    vVectorFields_(0),
    vSphericalTensorFields_(0),
    vSymmTensorFields_(0),
    vTensorFields_(0),

    sScalarFields_(0),
    sVectorFields_(0),
    sSphericalTensorFields_(0),
    sSymmTensorFields_(0),
    sTensorFields_(0)
{
    read(dict);
    forAll(fieldNames_, fieldi)
    {
        maxFieldNames_[fieldi] =
            IOobject::member(fieldNames_[fieldi])
          + "Max"
          + IOobject::group(fieldNames_[fieldi]);

        createMax<volScalarField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<volVectorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<volSphericalTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<volSymmTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<volTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );

        createMax<surfaceScalarField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<surfaceVectorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<surfaceSphericalTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<surfaceSymmTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
        createMax<surfaceTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi]
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMax::~fieldMax()
{
    clearOldFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);

    Log << endl;

    return true;
}


void Foam::functionObjects::fieldMax::updateMesh(const mapPolyMesh& mpm)
{
    forAll(fieldNames_, fieldi)
    {
        mapMax<volScalarField>
        (
            maxFieldNames_[fieldi], mpm, vScalarFields_
        );
        mapMax<volVectorField>
        (
            maxFieldNames_[fieldi], mpm, vVectorFields_
        );
        mapMax<volSphericalTensorField>
        (
            maxFieldNames_[fieldi], mpm, vSphericalTensorFields_
        );
        mapMax<volSymmTensorField>
        (
            maxFieldNames_[fieldi], mpm, vSymmTensorFields_
        );
        mapMax<volTensorField>
        (
            maxFieldNames_[fieldi], mpm, vTensorFields_
        );

        mapMax<surfaceScalarField>
        (
            maxFieldNames_[fieldi], mpm, sScalarFields_
        );
        mapMax<surfaceVectorField>
        (
            maxFieldNames_[fieldi], mpm, sVectorFields_
        );
        mapMax<surfaceSphericalTensorField>
        (
            maxFieldNames_[fieldi], mpm, sSphericalTensorFields_
        );
        mapMax<surfaceSymmTensorField>
        (
            maxFieldNames_[fieldi], mpm, sSymmTensorFields_
        );
        mapMax<surfaceTensorField>
        (
            maxFieldNames_[fieldi], mpm, sTensorFields_
        );
    }

    setOldFields(mpm);
}


bool Foam::functionObjects::fieldMax::execute()
{
    forAll(fieldNames_, fieldi)
    {
        updateMax<volScalarField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            vScalarFields_
        );
        updateMax<volVectorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            vVectorFields_
        );
        updateMax<volSphericalTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            vSphericalTensorFields_
        );
        updateMax<volSymmTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            vSymmTensorFields_
        );
        updateMax<volTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            vTensorFields_
        );

        updateMax<surfaceScalarField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            sScalarFields_
        );
        updateMax<surfaceVectorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            sVectorFields_
        );
        updateMax<surfaceSphericalTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            sSphericalTensorFields_
        );
        updateMax<surfaceSymmTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            sSymmTensorFields_
        );
        updateMax<surfaceTensorField>
        (
            fieldNames_[fieldi],
            maxFieldNames_[fieldi],
            sTensorFields_
        );
    }

    return true;
}


void Foam::functionObjects::fieldMax::clearOldFields()
{
    if (cellMap_)
    {
        delete cellMap_;
        cellMap_ = 0;
        delete rCellMap_;
        rCellMap_ = 0;

        // Clear toc and fields
        vScalarFields_.clear();
        vVectorFields_.clear();
        vSphericalTensorFields_.clear();
        vSymmTensorFields_.clear();
        vTensorFields_.clear();

        sScalarFields_.clear();
        sVectorFields_.clear();
        sSphericalTensorFields_.clear();
        sSymmTensorFields_.clear();
        sTensorFields_.clear();
    }
}


void Foam::functionObjects::fieldMax::setOldFields(const mapPolyMesh& mpm)
{
    clearOldFields();

    forAll(fieldNames_, fieldi)
    {
        createOld<volScalarField>
        (
            maxFieldNames_[fieldi],
            vScalarFields_
        );
        createOld<volVectorField>
        (
            maxFieldNames_[fieldi],
            vVectorFields_
        );
        createOld<volSphericalTensorField>
        (
            maxFieldNames_[fieldi],
            vSphericalTensorFields_
        );
        createOld<volSymmTensorField>
        (
            maxFieldNames_[fieldi],
            vSymmTensorFields_
        );
        createOld<volTensorField>
        (
            maxFieldNames_[fieldi],
            vTensorFields_
        );

        createOld<surfaceScalarField>
        (
            maxFieldNames_[fieldi],
            sScalarFields_
        );
        createOld<surfaceVectorField>
        (
            maxFieldNames_[fieldi],
            sVectorFields_
        );
        createOld<surfaceSphericalTensorField>
        (
            maxFieldNames_[fieldi],
            sSphericalTensorFields_
        );
        createOld<surfaceSymmTensorField>
        (
            maxFieldNames_[fieldi],
            sSymmTensorFields_
        );
        createOld<surfaceTensorField>
        (
            maxFieldNames_[fieldi],
            sTensorFields_
        );
    }

    cellMap_ = new labelList(mpm.cellMap());
    rCellMap_ = new labelList(mpm.reverseCellMap());
}


bool Foam::functionObjects::fieldMax::write()
{
    forAll(fieldNames_, fieldi)
    {
        writeField<volScalarField>(maxFieldNames_[fieldi]);
        writeField<volVectorField>(maxFieldNames_[fieldi]);
        writeField<volSphericalTensorField>(maxFieldNames_[fieldi]);
        writeField<volSymmTensorField>(maxFieldNames_[fieldi]);
        writeField<volTensorField>(maxFieldNames_[fieldi]);

        writeField<surfaceScalarField>(maxFieldNames_[fieldi]);
        writeField<surfaceVectorField>(maxFieldNames_[fieldi]);
        writeField<surfaceSphericalTensorField>(maxFieldNames_[fieldi]);
        writeField<surfaceSymmTensorField>(maxFieldNames_[fieldi]);
        writeField<surfaceTensorField>(maxFieldNames_[fieldi]);
    }
    return true;
}


// ************************************************************************* //
