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

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class InType, class OutType>
Foam::label Foam::Equation<InType, OutType>::nVar() const
{
    return nVar_;
}


template<class InType, class OutType>
Foam::label Foam::Equation<InType, OutType>::nEqns() const
{
    return nEqns_;
}


template<class InType, class OutType>
const InType& Foam::Equation<InType, OutType>::lower() const
{
    return lowerLimits_;
}


template<class InType, class OutType>
InType& Foam::Equation<InType, OutType>::lower()
{
    return lowerLimits_;
}


template<class InType, class OutType>
const InType& Foam::Equation<InType, OutType>::upper() const
{
    return upperLimits_;
}


template<class InType, class OutType>
InType& Foam::Equation<InType, OutType>::upper()
{
    return upperLimits_;
}


template<class InType, class OutType>
const InType& Foam::Equation<InType, OutType>::dx() const
{
    return dx_;
}


template<class InType, class OutType>
InType& Foam::Equation<InType, OutType>::dx()
{
    return dx_;
}


template<class InType, class OutType>
Foam::scalar Foam::Equation<InType, OutType>::dx(const label cmpti) const
{
    return component(dx_, cmpti);
}



template<class InType, class OutType>
void Foam::Equation<InType, OutType>::limit(InType& x) const
{
    x = max(lowerLimits_, min(upperLimits_, x));
}


// ************************************************************************* //
