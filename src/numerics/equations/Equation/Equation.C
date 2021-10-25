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

#include "Equation.H"

// * * * * * * * * * * * * * Private member functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class InType, class OutType>
Foam::Equation<InType, OutType>::Equation
(
    const label nVar,
    const label nEqns,
    const InType& lowerLimits,
    const InType& upperLimits
)
:
    nVar_(nVar),
    nEqns_(nEqns),
    lowerLimits_(lowerLimits),
    upperLimits_(upperLimits),
    dx_(lowerLimits)
{
    for (label cmpti = 0; cmpti < nVar; cmpti++)
    {
        setComponent(dx_, cmpti) = 1e-6;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class InType, class OutType>
Foam::Equation<InType, OutType>::~Equation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
