/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019-2020
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

#include "integrationSystem.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class fieldType>
void Foam::integrationSystem::storeOld
(
    const label stepi,
    fieldType& f,
    PtrList<fieldType>& fList
) const
{
    // Correct old field for mesh motion before storage
    if (f.mesh().moving() && stepi == 1)
    {
        f.ref() *= f.mesh().Vsc0()/f.mesh().Vsc();
    }

    // Store fields if needed later
    if (oldIs_[stepi - 1] != -1)
    {
        fList.set
        (
            oldIs_[stepi - 1],
            new fieldType(f.name() + Foam::name(stepi - 1), f)
        );
    }
}


template<class fieldType>
void Foam::integrationSystem::storeDelta
(
    const label stepi,
    const fieldType& f,
    PtrList<fieldType>& fList
) const
{
    // Store fields if needed later
    if (deltaIs_[stepi - 1] != -1)
    {
        fList.set
        (
            deltaIs_[stepi - 1],
            new fieldType(f.name() + Foam::name(stepi - 1), f)
        );
    }
}


template<class fieldType>
void Foam::integrationSystem::blendSteps
(
    const label stepi,
    const labelList& indices,
    fieldType& f,
    const PtrList<fieldType>& fList,
    const scalarList& scales
) const
{
    // Scale current step by weight
    f *= scales[stepi - 1];
    for (label i = 0; i < stepi - 1; i++)
    {
        label fi = indices[i];
        if (fi != -1 && scales[fi] != 0)
        {
            f += scales[fi]*fList[fi];
        }
    }
}


template<class fieldType>
void Foam::integrationSystem::blendOld
(
    const label stepi,
    fieldType& f,
    const PtrList<fieldType>& fList,
    const scalarList& scales
) const
{
    blendSteps(stepi, oldIs_, f, fList, scales);
}


template<class fieldType>
void Foam::integrationSystem::blendDelta
(
    const label stepi,
    fieldType& f,
    const PtrList<fieldType>& fList,
    const scalarList& scales
) const
{
    blendSteps(stepi, deltaIs_, f, fList, scales);
}


template<class fieldType>
void Foam::integrationSystem::clearOld(PtrList<fieldType>& fList) const
{
    fList.clear();
    fList.resize(nOld_);
}


template<class fieldType>
void Foam::integrationSystem::clearDelta(PtrList<fieldType>& fList) const
{
    fList.clear();
    fList.resize(nDelta_);
}
// ************************************************************************* //
