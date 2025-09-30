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

#include "regionSolver.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class FieldType>
Foam::scalar Foam::regionSolver::residual
(
    const FieldType& f,
    const FieldType& fOld
)
{
    scalar error = Zero;
    forAll(f, i)
    {
        error += magSqr(f[i] - fOld[i]);
    }
    reduce(error, sumOp<scalar>());
    return sqrt(error);
}


// ************************************************************************* //
