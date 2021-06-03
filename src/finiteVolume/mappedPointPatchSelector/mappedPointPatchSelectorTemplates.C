/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2020
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

template<class Type>
void Foam::mappedPointPatchSelector::distribute(List<Type>& lst) const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->distribute(lst);
    }

    NotImplemented;
    return mappedMovingPatchPtr_->distribute(lst);
}


template<class Type, class CombineOp>
void Foam::mappedPointPatchSelector::distribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->distribute(lst, cop);
    }

    NotImplemented;
    return mappedMovingPatchPtr_->distribute(lst, cop);
}


template<class Type>
void Foam::mappedPointPatchSelector::reverseDistribute(List<Type>& lst) const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->reverseDistribute(lst);
    }

    NotImplemented;
    return mappedMovingPatchPtr_->reverseDistribute(lst);
}


template<class Type, class CombineOp>
void Foam::mappedPointPatchSelector::reverseDistribute
(
    List<Type>& lst,
    const CombineOp& cop
) const
{
    if (mappedMovingPatchPtr_)
    {
        return mappedMovingPatchPtr_->reverseDistribute(lst, cop);
    }

    NotImplemented;
    return mappedMovingPatchPtr_->reverseDistribute(lst, cop);
}


// ************************************************************************* //
