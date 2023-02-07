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
#include "accelerationSchemeNew.H"

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
    PtrList<accelerationScheme>(patches.size()),
    fieldName_(fieldName),
    mesh_(mesh)
{
    Info<< "Selecting acceleration schemes for " << fieldName << endl;
    bool found = false;
    label nPatches = 0;
    #define setSchemes(Type, Patch, Mesh)                                  \
    if (mesh_.foundObject<GeometricField<Type, Patch, Mesh>>(fieldName))   \
    {                                                                      \
        found = true;                                                      \
        GeometricField<Type, Patch, Mesh>& field =                         \
            mesh_.lookupObjectRef<GeometricField<Type, Patch, Mesh>>       \
            (                                                              \
                fieldName                                                  \
            );                                                             \
        forAll(patches, pi)                                                \
        {                                                                  \
            const label patchi = patches[pi];                              \
            if (field.boundaryField()[patchi].fixesValue())                \
            {                                                              \
                this->set                                                  \
                (                                                          \
                    nPatches++,                                            \
                    accelerationScheme::New                                \
                    (                                                      \
                        field,                                             \
                        patchi,                                            \
                        dict                                               \
                    )                                                      \
                );                                                         \
            }                                                              \
            else                                                           \
            {                                                              \
                WarningInFunction                                          \
                    << "Trying to relax "                                  \
                    << field.mesh().boundary()[patchi].name()              \
                    << "for " << field.name() << " but it does not"        \
                    << "fix a value. Skipping" << endl;                    \
            }                                                              \
        }                                                                  \
        this->setSize(nPatches);                                           \
    }

    FOR_ALL_FIELD_TYPES(setSchemes, fvPatchField, volMesh);
    FOR_ALL_FIELD_TYPES(setSchemes, fvsPatchField, surfaceMesh);
    FOR_ALL_FIELD_TYPES(setSchemes, pointPatchField, pointMesh);

    #undef setSchemes

    if (!found)
    {
        FatalErrorInFunction
            << "Could not find " << fieldName << endl;
    }
    this->read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAccelerationScheme::~fieldAccelerationScheme()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fieldAccelerationScheme::readControls(const dictionary& dict)
{
    if (!dict.isDict(mesh_.name()) || !this->size())
    {
        return false;
    }
    const dictionary& regionDict =
        dict.subDict(mesh_.name());
    if (regionDict.isDict(fieldName_))
    {
        regionDict.subDict(fieldName_).lookup("tolerance") >> tolerance_;
        regionDict.subDict(fieldName_).lookup("relTol") >> relTol_;
        Info<< indent << "Tolerances for " << fieldName_ << " (abs/rel): "
            << tolerance_ << "/" << relTol_ << endl;
        return true;
    }
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


void Foam::fieldAccelerationScheme::clear()
{
    forAll(*this, i)
    {
        this->operator[](i).clear();
    }
}


Foam::scalar Foam::fieldAccelerationScheme::error() const
{
    scalar errorSqr = Zero;
    forAll(*this, i)
    {
        errorSqr = sqr(this->operator[](i).error());
    }
    return sqrt(errorSqr);
}


Foam::scalar Foam::fieldAccelerationScheme::relError() const
{
    scalar errorSqr = Zero;
    scalar initErrorSqr = Zero;
    forAll(*this, i)
    {
        errorSqr = sqr(this->operator[](i).error());
        initErrorSqr = sqr(this->operator[](i).initError());
    }
    return sqrt(errorSqr)/(sqrt(initErrorSqr) + small);
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
