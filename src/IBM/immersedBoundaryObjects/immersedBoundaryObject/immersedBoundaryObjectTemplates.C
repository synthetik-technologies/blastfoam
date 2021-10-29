/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021 Synthetik Applied Technology
     \\/     M anipulation  |
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

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::immersedBoundaryObject::interpolateTo(const Field<Type>& vf) const
{
    tmp<Field<Type>> tmpF(new Field<Type>(nFaces(), Zero));
    Field<Type>& f(tmpF.ref());
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(interpToCells_, i)
    {
        forAll(interpToCells_[i], j)
        {
            label celli = interpToCells_[i][j];
            if (celli >= 0)
            {
                f[i] += vf[celli]*interpToWeights_[i][j];
            }
        }
    }
    reduce(f, sumOp<List<Type>>());
    return tmpF;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::immersedBoundaryObject::patchInternalField(const Field<Type>& vf) const
{
    tmp<Field<Type>> tmpF(new Field<Type>(nFaces(), Zero));
    Field<Type>& f(tmpF.ref());
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(patchInternalCells_, i)
    {
        label celli = patchInternalCells_[i];
        if (celli >= 0)
        {
            f[i] = vf[celli];
        }
    }
    return tmpF;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::immersedBoundaryObject::patchExternalField(const Field<Type>& vf) const
{
    tmp<Field<Type>> tmpF(new Field<Type>(nFaces(), Zero));
    Field<Type>& f(tmpF.ref());
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(patchExternalCells_, i)
    {
        label celli = patchExternalCells_[i];
        if (celli >= 0)
        {
            f[i] = vf[celli];
        }
    }
    return tmpF;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::immersedBoundaryObject::boundaryValues(const Field<Type>& vf) const
{
    tmp<Field<Type>> tmpF(new Field<Type>(nFaces(), Zero));
    Field<Type>& f(tmpF.ref());
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*boundaryCellsPtr_, i)
    {
        f[i] = vf[(*boundaryCellsPtr_)[i]];
    }
    return tmpF;

}


template<class Type>
void Foam::immersedBoundaryObject::interpolateFrom
(
    const Field<Type>& sf,
    Field<Type>& vf,
    const bool clear
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    if (clear)
    {
        vf = Zero;
    }

    forAll(*shellCellsPtr_, i)
    {
        label celli = (*shellCellsPtr_)[i];
        forAll(interpFromPoints_[i], j)
        {
            vf[celli] +=
                sf[interpFromPoints_[i][j]]*interpFromWeights_[i][j];
        }
    }
}


template<class Type>
void Foam::immersedBoundaryObject::setInternal
(
    Field<Type>& vf,
    const Type& val
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*internalCellsPtr_, i)
    {
        vf[(*internalCellsPtr_)[i]] = val;
    }
}


template<class Type, class CombineOp>
void Foam::immersedBoundaryObject::setInternal
(
    Field<Type>& vf,
    const Type& val,
    const CombineOp& cop
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*internalCellsPtr_, i)
    {
        cop(vf[(*internalCellsPtr_)[i]], val);
    }
}


template<class Type>
void Foam::immersedBoundaryObject::setInternal
(
    Field<Type>& vf,
    const Field<Type>& vals
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*internalCellsPtr_, i)
    {
        vf[(*internalCellsPtr_)[i]] = vals[i];
    }
}


template<class Type, class CombineOp>
void Foam::immersedBoundaryObject::setInternal
(
    Field<Type>& vf,
    const Field<Type>& vals,
    const CombineOp& cop
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*internalCellsPtr_, i)
    {
        cop(vf[(*internalCellsPtr_)[i]], vals[i]);
    }
}


template<class Type>
void Foam::immersedBoundaryObject::setShell
(
    Field<Type>& vf,
    const Type& val
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*shellCellsPtr_, i)
    {
        vf[(*shellCellsPtr_)[i]] = val;
    }
}


template<class Type, class CombineOp>
void Foam::immersedBoundaryObject::setShell
(
    Field<Type>& vf,
    const Type& val,
    const CombineOp& cop
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*shellCellsPtr_, i)
    {
        cop(vf[(*shellCellsPtr_)[i]], val);
    }
}


template<class Type>
void Foam::immersedBoundaryObject::setBoundary
(
    Field<Type>& vf,
    const Type& val
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*boundaryCellsPtr_, i)
    {
        vf[(*boundaryCellsPtr_)[i]] = val;
    }
}


template<class Type, class CombineOp>
void Foam::immersedBoundaryObject::setBoundary
(
    Field<Type>& vf,
    const Type& val,
    const CombineOp& cop
) const
{
    if (!internalCellsPtr_)
    {
        calcMapping();
    }
    forAll(*boundaryCellsPtr_, i)
    {
        cop(vf[(*boundaryCellsPtr_)[i]], val);
    }
}


// ************************************************************************* //
