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
    for (label i = 0; i < nVar; i++)
    {
        setComponent(dx_, i) = 1e-6;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class InType, class OutType>
Foam::Equation<InType, OutType>::~Equation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class InType, class OutType>
Foam::tmp<Foam::scalarField>
Foam::Equation<InType, OutType>::lowerLimits() const 
{
    tmp<scalarField> tlower(new scalarField(nVar_));
    scalarField& lower = tlower.ref();
    for (label i = 0; i < nVar_; i++)
    {
        lower[i] = component(lowerLimits_, i);
    }
    return tlower;
}


template<class InType, class OutType>
Foam::tmp<Foam::scalarField>
Foam::Equation<InType, OutType>::upperLimits() const 
{
    tmp<scalarField> tupper(new scalarField(nVar_));
    scalarField& upper = tupper.ref();
    for (label i = 0; i < nVar_; i++)
    {
        upper[i] = component(upperLimits_, i);
    }
    return tupper;
}


template<class InType, class OutType>
Foam::tmp<Foam::scalarField>
Foam::Equation<InType, OutType>::dx() const 
{
    tmp<scalarField> tdx(new scalarField(nVar_));
    scalarField& dx = tdx.ref();
    for (label i = 0; i < nVar_; i++)
    {
        dx[i] = component(dx_, i);
    }
    return tdx;
}


// ************************************************************************* //
