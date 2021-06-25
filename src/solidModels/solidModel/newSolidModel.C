/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "solidModel.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidModel> Foam::solidModel::New(fvMesh& mesh)
{
    word solidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid model is created, otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidModel")
            >> solidModelTypeName;
    }

    Info<< "Selecting solidModel " << solidModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solidModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::New(Time&, const word&)"
        )   << "Unknown solidModel type " << solidModelTypeName
            << endl << endl
            << "Valid solidModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidModel>(cstrIter()(mesh));
}


Foam::autoPtr<Foam::solidModel> Foam::solidModel::NewLU(fvMesh& mesh)
{
    word solidModelTypeName;

    // Enclose the creation of the dictionary to ensure it is
    // deleted before the fluid model is created, otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary solidProperties
        (
            IOobject
            (
                "solidProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        solidProperties.lookup("solidModel")
            >> solidModelTypeName;
    }

    Info<< "Selecting solidModel " << solidModelTypeName << endl;

    lagrangianConstructorTable::iterator cstrIter =
        lagrangianConstructorTablePtr_->find(solidModelTypeName);

    if (cstrIter == lagrangianConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "solidModel::NewLU(Time&, const word&)"
        )   << "Unknown lagrangian solidModel type " << solidModelTypeName
            << endl << endl
            << "Valid lagrangian solidModel types are :" << endl
            << lagrangianConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<solidModel>(cstrIter()(mesh));
}
