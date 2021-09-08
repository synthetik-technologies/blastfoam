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

#include "fieldMinMax.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldMinMax, 0);
    addToRunTimeSelectionTable(functionObject, fieldMinMax, dictionary);

    // Add old name
    addNamedToRunTimeSelectionTable
    (
        functionObject,
        fieldMinMax,
        dictionary,
        fieldMax
    );
}
}

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    3
>::names[] ={"component", "magnitude", "componentMag"};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::minMaxType,
    2
>::names[] = {"min", "max"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::modeType,
    3
> Foam::functionObjects::fieldMinMax::modeTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMinMax::minMaxType,
    2
> Foam::functionObjects::fieldMinMax::minMaxTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::fieldMinMax
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    restartOnRestart_(dict.lookupOrDefault("restartOnRestart", false)),
    mode_(modeType::cmpt),
    minMax_(minMaxType::max),
    fieldNames_(dict.lookup("fields")),
    computeFieldNames_(fieldNames_.size()),

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

    word minMaxName;
    if (minMax_ == minMaxType::min)
    {
        minMaxName = "Min";
    }
    else
    {
        minMaxName = "Max";
    }

    if (mode_ == modeType::mag)
    {
        minMaxName += "Mag";
    }
    else if (mode_ == modeType::cmptMag)
    {
        minMaxName += "CmptMag";
    }

    forAll(fieldNames_, fieldi)
    {
        computeFieldNames_[fieldi] =
            IOobject::member(fieldNames_[fieldi])
          + minMaxName
          + IOobject::group(fieldNames_[fieldi]);

        createMinMax<volScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<volVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<volSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<volSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<volTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );

        createMinMax<surfaceScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<surfaceVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<surfaceSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<surfaceSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMinMax<surfaceTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldMinMax::~fieldMinMax()
{
    clearOldFields();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldMinMax::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    Log << type() << " " << name() << ":" << nl;

    dict.readIfPresent("restartOnRestart", restartOnRestart_);
    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "component")];
    minMax_ = minMaxTypeNames_[dict.lookupOrDefault<word>("minMax", "max")];

    Log << endl;

    return true;
}


void Foam::functionObjects::fieldMinMax::updateMesh(const mapPolyMesh& mpm)
{
    forAll(fieldNames_, fieldi)
    {
        map<volScalarField>
        (
            computeFieldNames_[fieldi], mpm, vScalarFields_
        );
        map<volVectorField>
        (
            computeFieldNames_[fieldi], mpm, vVectorFields_
        );
        map<volSphericalTensorField>
        (
            computeFieldNames_[fieldi], mpm, vSphericalTensorFields_
        );
        map<volSymmTensorField>
        (
            computeFieldNames_[fieldi], mpm, vSymmTensorFields_
        );
        map<volTensorField>
        (
            computeFieldNames_[fieldi], mpm, vTensorFields_
        );

        map<surfaceScalarField>
        (
            computeFieldNames_[fieldi], mpm, sScalarFields_
        );
        map<surfaceVectorField>
        (
            computeFieldNames_[fieldi], mpm, sVectorFields_
        );
        map<surfaceSphericalTensorField>
        (
            computeFieldNames_[fieldi], mpm, sSphericalTensorFields_
        );
        map<surfaceSymmTensorField>
        (
            computeFieldNames_[fieldi], mpm, sSymmTensorFields_
        );
        map<surfaceTensorField>
        (
            computeFieldNames_[fieldi], mpm, sTensorFields_
        );
    }

    setOldFields(mpm);
}


bool Foam::functionObjects::fieldMinMax::execute()
{
    forAll(fieldNames_, fieldi)
    {
        update<volScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vScalarFields_
        );
        update<volVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vVectorFields_
        );
        update<volSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vSphericalTensorFields_
        );
        update<volSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vSymmTensorFields_
        );
        update<volTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vTensorFields_
        );

        update<surfaceScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sScalarFields_
        );
        update<surfaceVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sVectorFields_
        );
        update<surfaceSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sSphericalTensorFields_
        );
        update<surfaceSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sSymmTensorFields_
        );
        update<surfaceTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sTensorFields_
        );
    }

    return true;
}


void Foam::functionObjects::fieldMinMax::clearOldFields()
{
    if (cellMap_.valid())
    {
        cellMap_.clear();
        rCellMap_.clear();

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


void Foam::functionObjects::fieldMinMax::setOldFields(const mapPolyMesh& mpm)
{
    clearOldFields();

    forAll(fieldNames_, fieldi)
    {
        createOld<volScalarField>
        (
            computeFieldNames_[fieldi],
            vScalarFields_
        );
        createOld<volVectorField>
        (
            computeFieldNames_[fieldi],
            vVectorFields_
        );
        createOld<volSphericalTensorField>
        (
            computeFieldNames_[fieldi],
            vSphericalTensorFields_
        );
        createOld<volSymmTensorField>
        (
            computeFieldNames_[fieldi],
            vSymmTensorFields_
        );
        createOld<volTensorField>
        (
            computeFieldNames_[fieldi],
            vTensorFields_
        );

        createOld<surfaceScalarField>
        (
            computeFieldNames_[fieldi],
            sScalarFields_
        );
        createOld<surfaceVectorField>
        (
            computeFieldNames_[fieldi],
            sVectorFields_
        );
        createOld<surfaceSphericalTensorField>
        (
            computeFieldNames_[fieldi],
            sSphericalTensorFields_
        );
        createOld<surfaceSymmTensorField>
        (
            computeFieldNames_[fieldi],
            sSymmTensorFields_
        );
        createOld<surfaceTensorField>
        (
            computeFieldNames_[fieldi],
            sTensorFields_
        );
    }

    cellMap_.set(new labelList(mpm.cellMap()));
    rCellMap_.set(new labelList(mpm.reverseCellMap()));
}


bool Foam::functionObjects::fieldMinMax::write()
{
    if (obr_.time().timeIndex() == 0)
    {
        return false;
    }

    forAll(fieldNames_, fieldi)
    {
        writeField<volScalarField>(computeFieldNames_[fieldi]);
        writeField<volVectorField>(computeFieldNames_[fieldi]);
        writeField<volSphericalTensorField>(computeFieldNames_[fieldi]);
        writeField<volSymmTensorField>(computeFieldNames_[fieldi]);
        writeField<volTensorField>(computeFieldNames_[fieldi]);

        writeField<surfaceScalarField>(computeFieldNames_[fieldi]);
        writeField<surfaceVectorField>(computeFieldNames_[fieldi]);
        writeField<surfaceSphericalTensorField>(computeFieldNames_[fieldi]);
        writeField<surfaceSymmTensorField>(computeFieldNames_[fieldi]);
        writeField<surfaceTensorField>(computeFieldNames_[fieldi]);
    }
    return true;
}


// ************************************************************************* //
