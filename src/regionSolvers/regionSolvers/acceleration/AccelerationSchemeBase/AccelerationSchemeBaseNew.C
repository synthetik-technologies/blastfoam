/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Synthetik Applied Technologies
     \\/     M anipulation  |
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
------------------------------------------------------------------------*/

#include "AccelerationSchemeBase.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

template<class Type, template<class> class Patch, class Mesh>
Foam::autoPtr<Foam::accelerationScheme>
Foam::AccelerationSchemeBase<Type, Patch, Mesh>::New
(
    GeometricField<Type, Patch, Mesh>& field,
    const label patchi,
    const dictionary& dict
)
{
    incrIndent(Info);
    Info<< indent << field.mesh().boundary()[patchi].name() << ": ";
    word accelerationType =
        accelerationScheme::schemesDict
        (
            dict,
            field().mesh().thisDb().name(),
            field.name(),
            field().mesh().boundary()[patchi].name()
        ).template lookupOrDefault<word>(accelerationScheme::typeName, "none");

    Info<< accelerationType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(accelerationType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown acceleration scheme type "
            << accelerationType << endl << endl
            << "Valid acceleration schemes are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    autoPtr<Foam::accelerationScheme> scheme(cstrIter()(field, patchi, dict));
    Info<< endl << decrIndent;
    return scheme;
}


// ************************************************************************* //
