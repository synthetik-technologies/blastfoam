/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2020-04-02 Jeff Heylmun:    Modified class for a density based thermodynamic
                            class
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

\*---------------------------------------------------------------------------*/

#include "solidBlastThermo.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "coordinateSystem.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidBlastThermo, 0);
    defineRunTimeSelectionTable(solidBlastThermo, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidBlastThermo::solidBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    physicalProperties(mesh, phaseName),
    blastThermo(mesh, dict, phaseName)
{
    this->properties().dictionary::operator=(dict);
}


void Foam::solidBlastThermo::initializeFields()
{
    if (!this->isotropic())
    {
        KappaPtr_.set
        (
            new volVectorField
            (
                IOobject
                (
                    basicThermo::phasePropertyName("Kappa"),
                    mesh().time().name(),
                    mesh()
                ),
                mesh(),
                dimensionedVector(kappa_.dimensions(), Zero)
            )
        );
    }
    updateRho();
    if (!e_.headerOk())
    {
        e_ == this->calce();
    }

    correct();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidBlastThermo> Foam::solidBlastThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    return blastThermo::New<solidBlastThermo>
    (
        mesh,
        dict.optionalSubDict("mixture"),
        phaseName,
        phaseName
    );
}


Foam::autoPtr<Foam::solidBlastThermo> Foam::solidBlastThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    const IOdictionary dict
    (
        physicalProperties::findModelDict(mesh, phaseName)
    );

    return blastThermo::New<solidBlastThermo>
    (
        mesh,
        dict.optionalSubDict("mixture"),
        phaseName,
        phaseName
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidBlastThermo::~solidBlastThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidBlastThermo::nu() const
{
    return volScalarField::New
    (
        "nu",
        this->T_.mesh(),
        dimKinematicViscosity
    );
}

// ************************************************************************* //
