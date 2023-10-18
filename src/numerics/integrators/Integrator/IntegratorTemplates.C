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

#include "Integrator.H"
#include "Simpson13Integrator.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type, class Adapt>
Foam::autoPtr<Foam::Integrator<Type, Adapt>> Foam::Integrator<Type, Adapt>::New
(
    const equationType& eqn,
    const dictionary& dict
)
{
    return New
    (
        dict.lookupOrDefault<word>
        (
            "integrator",
            Simpson13Integrator<Type, Adapt>::typeName
        ),
        eqn,
        dict
    );
}


template<class Type, class Adapt>
Foam::autoPtr<Foam::Integrator<Type, Adapt>> Foam::Integrator<Type, Adapt>::New
(
    const word& integratorTypeName,
    const equationType& eqn,
    const dictionary& dict
)
{
    Info<< "Selecting integrator: " << integratorTypeName << endl;
    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(integratorTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown integrator type "
            << integratorTypeName << nl << nl
            << "Valid integrators are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Integrator<Type, Adapt>>(cstrIter()(eqn, dict));
}


template<class Type, class Adapt>
Foam::autoPtr<Foam::Integrator<Type, Adapt>> Foam::Integrator<Type, Adapt>::New
(
    const equationType& eqn,
    const integrator& inter
)
{
    typename inputsConstructorTable::iterator cstrIter =
        inputsConstructorTablePtr_->find(inter.type());

    if (cstrIter == inputsConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown integrator type "
            << inter.type() << nl << nl
            << "Valid integrators are : " << endl
            << inputsConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<Integrator<Type, Adapt>>(cstrIter()(eqn, inter));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class Adapt>
Foam::Integrator<Type, Adapt>::Integrator
(
    const equationType& eqn,
    const dictionary& dict
)
:
    integrator(dict),
    eqnPtr_(&eqn)
{}


template<class Type, class Adapt>
Foam::Integrator<Type, Adapt>::Integrator
(
    const equationType& eqn,
    const integrator& inter
)
:
    integrator(inter),
    eqnPtr_(&eqn)
{}

// ************************************************************************* //
