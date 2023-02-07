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

#include "accelerationSchemeList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(accelerationSchemeList, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::accelerationSchemeList::accelerationSchemeList
(
    const word& fieldName,
    const fvMesh& mesh,
    const labelList& patches,
    const dictionary& dict
)
:
    PtrList<accelerationScheme>(patches.size()),
    mesh_(mesh)
{
    Info<< "Selecting acceleration schemes for " << fieldName << endl;
    bool found = false;
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
            this->set                                                      \
            (                                                              \
                pi,                                                        \
                accelerationScheme::New                                    \
                (                                                          \
                    field,                                                 \
                    patchi,                                                \
                    dict                                                   \
                )                                                          \
            );                                                             \
        }                                                                  \
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
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::accelerationSchemeList::~accelerationSchemeList()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::accelerationSchemeList::updateError()
{
    forAll(*this, i)
    {
        this->operator[](i).updateError();
    }
}


void Foam::accelerationSchemeList::relax(const label iter)
{
    forAll(*this, i)
    {
        this->operator[](i).relax(iter);
    }
}


void Foam::accelerationSchemeList::clear()
{
    forAll(*this, i)
    {
        this->operator[](i).clear();
    }
}


Foam::scalar Foam::accelerationSchemeList::error() const
{
    scalar errorSqr = Zero;
    forAll(*this, i)
    {
        errorSqr = sqr(this->operator[](i).error());
    }
    return sqrt(errorSqr);
}


Foam::scalar Foam::accelerationSchemeList::relError() const
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


void Foam::accelerationSchemeList::read(const dictionary& dict)
{
    forAll(*this, i)
    {
        this->operator[](i).read(dict);
    }
}



// ************************************************************************* //
