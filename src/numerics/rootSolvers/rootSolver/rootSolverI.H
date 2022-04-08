/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2021
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

inline const Foam::scalarField& Foam::rootSolver::xTols() const
{
    return xTols_;
}


inline Foam::scalarField& Foam::rootSolver::xTols()
{
    return xTols_;
}


inline const Foam::scalarField& Foam::rootSolver::xAbsTols() const
{
    return xAbsTols_;
}


inline Foam::scalarField& Foam::rootSolver::xAbsTols()
{
    return xAbsTols_;
}


inline const Foam::scalarField& Foam::rootSolver::yTols() const
{
    return yTols_;
}


inline Foam::scalarField& Foam::rootSolver::yTols()
{
    return yTols_;
}


inline Foam::label Foam::rootSolver::nSteps() const
{
    return stepi_;
}


inline const Foam::scalarList& Foam::rootSolver::xErrors() const
{
    return xErrors_;
}


inline const Foam::scalarList& Foam::rootSolver::yErrors() const
{
    return yErrors_;
}


// ************************************************************************* //
