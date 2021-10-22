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

inline Foam::scalarList& Foam::multivariateMinimizationScheme::tolerances()
{
    return tolerances_;
}


inline Foam::label Foam::multivariateMinimizationScheme::nSteps() const
{
    return stepi_;
}


inline Foam::scalarList& Foam::multivariateMinimizationScheme::errors() const
{
    return errors_;
}


inline Foam::scalar Foam::multivariateMinimizationScheme::error() const
{
    return error_;
}

// ************************************************************************* //
