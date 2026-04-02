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

#include "timeIntegrationSystem.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FieldType>
void Foam::timeIntegrationSystem::storeOld
(
    FieldType& f,
    const bool conservative
)
{
    if (conservative)
    {
        fvTimeInt_->conservativeFieldsRef().insert(f.name());
    }
    if (storeOld_ && timeInt_->firstStep() && !timeInt_->restart())
    {
        f.storeOldTimes();
    }

    storeOld
    (
        f,
        fvTimeInt_->oldFieldsRef(f),
        conservative
    );
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeDelta(const FieldType& f)
{
    storeDelta(f(), fvTimeInt_->deltaFieldsRef(f));
}


template<class FieldType>
void Foam::timeIntegrationSystem::blendOld
(
    FieldType& f,
    const bool conservative
) const
{
    blendOld
    (
        f,
        fvTimeInt_->oldFields(f),
        conservative
    );
}


template<class FieldType>
void Foam::timeIntegrationSystem::blendDelta(FieldType& f) const
{
    blendDelta(f, fvTimeInt_->deltaFields(f));
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeAndBlendOld
(
    FieldType& f,
    const bool conservative
)
{
    if (conservative)
    {
        fvTimeInt_->conservativeFieldsRef().insert(f.name());
    }
    storeAndBlendOld
    (
        f,
        fvTimeInt_->oldFieldsRef(f),
        conservative
    );
}


template<class FieldType>
void Foam::timeIntegrationSystem::storeAndBlendDelta(FieldType& f)
{
    storeAndBlendDelta
    (
        f,
        fvTimeInt_->deltaFieldsRef(f)
    );
}


template<class FieldType>
Foam::tmp<FieldType> Foam::timeIntegrationSystem::calcDelta
(
    const FieldType& f
) const
{
    return calcDelta(f, fvTimeInt_->deltaFields(f));
}


template<class FieldType>
Foam::tmp<FieldType> Foam::timeIntegrationSystem::calcAndStoreDelta
(
    const FieldType& f
)
{
   return calcAndStoreDelta(f, fvTimeInt_->deltaFieldsRef(f));
}


template<class FieldType>
void Foam::timeIntegrationSystem::addOldField(const FieldType& f)
{
    fvTimeInt_->addOldField(f);
}


template<class FieldType>
void Foam::timeIntegrationSystem::addDeltaField(const FieldType& f)
{
    fvTimeInt_->addDeltaField(f);
}


template<class FieldType>
void Foam::timeIntegrationSystem::addOldDeltaField(const FieldType& f)
{
    fvTimeInt_->addOldField(f);
    fvTimeInt_->addDeltaField(f);
}

// ************************************************************************* //
