/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021-2022
     \\/     M anipulation  | Synthetik Applied Technologies
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

inline Foam::scalar Foam::univariateRootSolver::xError() const
{
    return xErrors_[0];
}


inline Foam::scalar Foam::univariateRootSolver::yError() const
{
    return yErrors_[0];
}


inline Foam::scalar Foam::univariateRootSolver::xTol() const
{
    return xTols_[0];
}


inline Foam::scalar& Foam::univariateRootSolver::xTol()
{
    return xTols_[0];
}


inline Foam::scalar Foam::univariateRootSolver::xAbsTol() const
{
    return xAbsTols_[0];
}


inline Foam::scalar& Foam::univariateRootSolver::xAbsTol()
{
    return xAbsTols_[0];
}


inline Foam::scalar Foam::univariateRootSolver::yTol() const
{
    return yTols_[0];
}


inline Foam::scalar& Foam::univariateRootSolver::yTol()
{
    return yTols_[0];
}


// ************************************************************************* //
