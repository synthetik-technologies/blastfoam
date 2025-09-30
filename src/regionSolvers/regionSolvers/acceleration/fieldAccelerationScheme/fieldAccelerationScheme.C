/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2021
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

#include "fieldAccelerationScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldAccelerationScheme, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAccelerationScheme::fieldAccelerationScheme
(
    const word& fieldName,
    const fvMesh& mesh,
    const labelList& patches,
    const dictionary& dict
)
:
    PtrList<accelerationScheme>(0),
    fieldName_(fieldName),
    mesh_(mesh),
    tolerance_(-1.0),
    relTol_(0.0)
{
    Info<< "Selecting acceleration schemes for " << fieldName << endl;
    incrIndent(Info);
    bool found = false;

    #define setAccelerators(Type, Patch, Mesh)                             \
    found =                                                                \
        setPatchAccelerators<Type, Patch, Mesh>                            \
        (                                                                  \
            patches,                                                       \
            dict                                                           \
        );                                                                 \

    FOR_ALL_FIELD_TYPES(setAccelerators, fvPatchField, volMesh);
    FOR_ALL_FIELD_TYPES(setAccelerators, fvsPatchField, surfaceMesh);
    FOR_ALL_FIELD_TYPES(setAccelerators, pointPatchField, pointMesh);

    #undef setAccelerators

    if (!found)
    {
        FatalErrorInFunction
            << "Could not find " << fieldName << endl;
    }
    this->read(dict);

    Info<< decrIndent << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAccelerationScheme::~fieldAccelerationScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fieldAccelerationScheme::readControls(const dictionary& dict)
{
    if
    (
        (
            mesh_.name() != polyMesh::defaultRegion
         && !dict.isDict(mesh_.name())
        )
     || !this->size())
    {
        tolerance_ = -1.0;
        return false;
    }
    const dictionary& regionDict =
        mesh_.name() == polyMesh::defaultRegion
      ? dict
      : dict.subDict(mesh_.name());
    if (!regionDict.isDict("outerCorrectorResidualControl"))
    {
        tolerance_ = -1.0;
        return false;
    }

    const dictionary& residualDict =
        regionDict.subDict("outerCorrectorResidualControl");

    if (residualDict.isDict(fieldName_))
    {
        residualDict.subDict(fieldName_).lookup("tolerance") >> tolerance_;
        residualDict.subDict(fieldName_).lookup("relTol") >> relTol_;
        Info<< indent << "Tolerances for " << fieldName_ << " (abs/rel): "
            << tolerance_ << "/" << relTol_ << endl;
        return true;
    }
    tolerance_ = -1.0;
    return false;
}


void Foam::fieldAccelerationScheme::updateError()
{
    forAll(*this, i)
    {
        this->operator[](i).updateError();
    }
}


void Foam::fieldAccelerationScheme::relax(const label iter)
{
    forAll(*this, i)
    {
        this->operator[](i).relax(iter);
    }
}


void Foam::fieldAccelerationScheme::clear(const bool full)
{
    forAll(*this, i)
    {
        this->operator[](i).clear(full);
    }
}


Foam::scalar Foam::fieldAccelerationScheme::error() const
{
    scalar maxError = Zero;
    forAll(*this, i)
    {
        maxError = max(maxError, this->operator[](i).error());
    }
    return maxError;
}


Foam::scalar Foam::fieldAccelerationScheme::relError() const
{
    scalar maxRelError = Zero;
    forAll(*this, i)
    {
        maxRelError = max(maxRelError, this->operator[](i).relError());
    }
    return maxRelError;
}


void Foam::fieldAccelerationScheme::storePrevIter()
{
    #define storePrevIterType(Type, Patch, Mesh)                           \
    if (storePrevIter<Type, Patch, Mesh>())                                \
    {                                                                      \
        return;                                                            \
    }

    FOR_ALL_FIELD_TYPES(storePrevIterType, fvPatchField, volMesh);
    FOR_ALL_FIELD_TYPES(storePrevIterType, fvsPatchField, surfaceMesh);
    FOR_ALL_FIELD_TYPES(storePrevIterType, pointPatchField, pointMesh);

    #undef storePrevIterType
}


void Foam::fieldAccelerationScheme::read(const dictionary& dict)
{
    forAll(*this, i)
    {
        this->operator[](i).read(dict);
    }
    readControls(dict);
}


void Foam::fieldAccelerationScheme::print(Ostream& os) const
{
    if (this->size())
    {
        os  <<fieldName_ << " error (abs/rel): "
            << error() << "/" << relError() << endl;
    }
}


// ************************************************************************* //
