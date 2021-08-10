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

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMax::modeType,
    3
>::names[] ={"component", "magnitude", "componentMag"};

template<>
const char* Foam::NamedEnum
<
    Foam::functionObjects::fieldMax::minMaxType,
    2
>::names[] = {"min", "max"};

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMax::modeType,
    3
> Foam::functionObjects::fieldMax::modeTypeNames_;

const Foam::NamedEnum
<
    Foam::functionObjects::fieldMax::minMaxType,
    2
> Foam::functionObjects::fieldMax::minMaxTypeNames_;


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

        createMax<volScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<volVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<volSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<volSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<volTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );

        createMax<surfaceScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<surfaceVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<surfaceSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<surfaceSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
        );
        createMax<surfaceTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi]
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
    mode_ = modeTypeNames_[dict.lookupOrDefault<word>("mode", "component")];
    minMax_ = minMaxTypeNames_[dict.lookupOrDefault<word>("minMax", "max")];

    Log << endl;

    return true;
}


void Foam::functionObjects::fieldMax::updateMesh(const mapPolyMesh& mpm)
{
    forAll(fieldNames_, fieldi)
    {
        mapMax<volScalarField>
        (
            computeFieldNames_[fieldi], mpm, vScalarFields_
        );
        mapMax<volVectorField>
        (
            computeFieldNames_[fieldi], mpm, vVectorFields_
        );
        mapMax<volSphericalTensorField>
        (
            computeFieldNames_[fieldi], mpm, vSphericalTensorFields_
        );
        mapMax<volSymmTensorField>
        (
            computeFieldNames_[fieldi], mpm, vSymmTensorFields_
        );
        mapMax<volTensorField>
        (
            computeFieldNames_[fieldi], mpm, vTensorFields_
        );

        mapMax<surfaceScalarField>
        (
            computeFieldNames_[fieldi], mpm, sScalarFields_
        );
        mapMax<surfaceVectorField>
        (
            computeFieldNames_[fieldi], mpm, sVectorFields_
        );
        mapMax<surfaceSphericalTensorField>
        (
            computeFieldNames_[fieldi], mpm, sSphericalTensorFields_
        );
        mapMax<surfaceSymmTensorField>
        (
            computeFieldNames_[fieldi], mpm, sSymmTensorFields_
        );
        mapMax<surfaceTensorField>
        (
            computeFieldNames_[fieldi], mpm, sTensorFields_
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
            computeFieldNames_[fieldi],
            vScalarFields_
        );
        updateMax<volVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vVectorFields_
        );
        updateMax<volSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vSphericalTensorFields_
        );
        updateMax<volSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vSymmTensorFields_
        );
        updateMax<volTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            vTensorFields_
        );

        updateMax<surfaceScalarField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sScalarFields_
        );
        updateMax<surfaceVectorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sVectorFields_
        );
        updateMax<surfaceSphericalTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sSphericalTensorFields_
        );
        updateMax<surfaceSymmTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sSymmTensorFields_
        );
        updateMax<surfaceTensorField>
        (
            fieldNames_[fieldi],
            computeFieldNames_[fieldi],
            sTensorFields_
        );
    }

    return true;
}


void Foam::functionObjects::fieldMax::clearOldFields()
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


void Foam::functionObjects::fieldMax::setOldFields(const mapPolyMesh& mpm)
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


bool Foam::functionObjects::fieldMax::write()
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
