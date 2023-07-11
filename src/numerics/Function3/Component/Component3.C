/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
03-12-2021 Synthetik Applied Technologies : Added Function3
-------------------------------------------------------------------------------
License
    This file is a derived work of OpenFOAM.

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

#include "Component3.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::Function3s::Component<Type>::read(const dictionary& dict)
{
    xFunc_ =
        dict.found("x")
      ? Function1<Type>::New("x", dict)
      : autoPtr<Function1<Type>>
        (
            new Function1s::Constant<Type>("x", Zero)
        );
    yFunc_ =
        dict.found("y")
      ? Function1<Type>::New("y", dict)
      : autoPtr<Function1<Type>>
        (
            new Function1s::Constant<Type>("y", Zero)
        );
    zFunc_ =
        dict.found("z")
      ? Function1<Type>::New("z", dict)
      : autoPtr<Function1<Type>>
        (
            new Function1s::Constant<Type>("z", Zero)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function3s::Component<Type>::Component
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction3<Type, Component<Type>>(name)
{
    read(dict);
}


template<class Type>
Foam::Function3s::Component<Type>::Component(const Component<Type>& se)
:
    FieldFunction3<Type, Component<Type>>(se),
    xFunc_(se.xFunc_, false),
    yFunc_(se.yFunc_, false),
    zFunc_(se.zFunc_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function3s::Component<Type>::~Component()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function3s::Component<Type>::write(Ostream& os) const
{
    writeEntry(os, xFunc_());
    writeEntry(os, yFunc_());
    writeEntry(os, zFunc_());
}


// ************************************************************************* //
