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

#include "regEquation.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class BaseEquation>
Foam::regEquation<Type, BaseEquation>::regEquation
(
    const objectRegistry& obr,
    const dictionary& dict
)
:
    regIOobject
    (
        IOobject
        (
            dict.lookup<word>("name"),
            obr.instance(),
            obr
        )
    ),
    BaseEquation<Type>(dict),
    obr_(obr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type, template<class> class BaseEquation>
Foam::regEquation<Type, BaseEquation>::~regEquation()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type, template<class> class BaseEquation>
Foam::autoPtr<typename BaseEquation<Type>::BaseEquation>
Foam::regEquation<Type, BaseEquation>::New
(
    const word& type,
    const objectRegistry& obr,
    const dictionary& dict
)
{
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type "
            << type <<  nl << nl
            << "Valid " << typeName << " types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalError);
    }

    autoPtr<regEquation<Type, BaseEquation>> eqnPtr(cstrIter()(obr, dict));
    return autoPtr<typename BaseEquation<Type>::BaseEquation>(eqnPtr.ptr());
}

// ************************************************************************* //
