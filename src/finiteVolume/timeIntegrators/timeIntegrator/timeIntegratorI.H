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

namespace Foam
{

// * * * * * * * * * * * * * * * Public Functions  * * * * * * * * * * * * * //

// volScalarFields
void timeIntegrator::addOldField(const volScalarField& f) const
{
    insertOldList(oldScalarFields_, f);
}


void timeIntegrator::addDeltaField(const volScalarField& f) const
{
    insertDeltaList(oldScalarFields_, f);
}


const PtrList<volScalarField>& timeIntegrator::oldFields
(
    const volScalarField& f
) const
{
    return lookupOld(oldScalarFields_, f);
}


PtrList<volScalarField>& timeIntegrator::oldFields
(
    const volScalarField& f
)
{
    return lookupOld(oldScalarFields_, f);
}


const PtrList<volScalarField>& timeIntegrator::deltaFields
(
    const volScalarField& f
) const
{
    return lookupDelta(deltaScalarFields_, f);
}


PtrList<volScalarField>& timeIntegrator::deltaFields
(
    const volScalarField& f
)
{
    return lookupDelta(deltaScalarFields_, f);
}


// volVectorField
void timeIntegrator::addOldField(const volVectorField& f) const
{
    insertOldList(oldVectorFields_, f);
}


void timeIntegrator::addDeltaField(const volVectorField& f) const
{
    insertDeltaList(deltaVectorFields_, f);
}


const PtrList<volVectorField>& timeIntegrator::oldFields
(
    const volVectorField& f
) const
{
    return lookupOld(oldVectorFields_, f);
}


PtrList<volVectorField>& timeIntegrator::oldFields
(
    const volVectorField& f
)
{
    return lookupOld(oldVectorFields_, f);
}

const PtrList<volVectorField>& timeIntegrator::deltaFields
(
    const volVectorField& f
) const
{
    return lookupDelta(deltaVectorFields_, f);
}


PtrList<volVectorField>& timeIntegrator::deltaFields
(
    const volVectorField& f
)
{
    return lookupDelta(deltaVectorFields_, f);
}


// volSphericalTensorField
void timeIntegrator::addOldField(const volSphericalTensorField& f) const
{
    insertOldList(oldSphTensorFields_, f);
}


void timeIntegrator::addDeltaField(const volSphericalTensorField& f) const
{
    insertDeltaList(deltaSphTensorFields_, f);
}


const PtrList<volSphericalTensorField>& timeIntegrator::oldFields
(
    const volSphericalTensorField& f
) const
{
    return lookupOld(oldSphTensorFields_, f);
}


PtrList<volSphericalTensorField>& timeIntegrator::oldFields
(
    const volSphericalTensorField& f
)
{
    return lookupOld(oldSphTensorFields_, f);
}

const PtrList<volSphericalTensorField>& timeIntegrator::deltaFields
(
    const volSphericalTensorField& f
) const
{
    return lookupDelta(deltaSphTensorFields_, f);
}


PtrList<volSphericalTensorField>& timeIntegrator::deltaFields
(
    const volSphericalTensorField& f
)
{
    return lookupDelta(deltaSphTensorFields_, f);
}


// volSymmTensorField
void timeIntegrator::addOldField(const volSymmTensorField& f) const
{
    insertOldList(oldSymmTensorFields_, f);
}


void timeIntegrator::addDeltaField(const volSymmTensorField& f) const
{
    insertDeltaList(deltaSymmTensorFields_, f);
}


const PtrList<volSymmTensorField>& timeIntegrator::oldFields
(
    const volSymmTensorField& f
) const
{
    return lookupOld(oldSymmTensorFields_, f);
}


PtrList<volSymmTensorField>& timeIntegrator::oldFields
(
    const volSymmTensorField& f
)
{
    return lookupOld(oldSymmTensorFields_, f);
}


const PtrList<volSymmTensorField>& timeIntegrator::deltaFields
(
    const volSymmTensorField& f
) const
{
    return lookupDelta(deltaSymmTensorFields_, f);
}


PtrList<volSymmTensorField>& timeIntegrator::deltaFields
(
    const volSymmTensorField& f
)
{
    return lookupDelta(deltaSymmTensorFields_, f);
}


// volTensorField
void timeIntegrator::addOldField(const volTensorField& f) const
{
    insertOldList(oldTensorFields_, f);
}


void timeIntegrator::addDeltaField(const volTensorField& f) const
{
    insertDeltaList(deltaTensorFields_, f);
}


const PtrList<volTensorField>& timeIntegrator::oldFields
(
    const volTensorField& f
) const
{
    return lookupOld(oldTensorFields_, f);
}


PtrList<volTensorField>& timeIntegrator::oldFields
(
    const volTensorField& f
)
{
    return lookupOld(oldTensorFields_, f);
}

const PtrList<volTensorField>& timeIntegrator::deltaFields
(
    const volTensorField& f
) const
{
    return lookupDelta(deltaTensorFields_, f);
}


PtrList<volTensorField>& timeIntegrator::deltaFields
(
    const volTensorField& f
)
{
    return lookupDelta(deltaTensorFields_, f);
}

} // End namespace Foam

// ************************************************************************* //
