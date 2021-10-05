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

#include "phaseBlastThermo.H"
#include "blastThermo.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(phaseBlastThermo, 0);
    defineRunTimeSelectionTable(phaseBlastThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phaseBlastThermo::phaseBlastThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
:
    timeIntegrationSystem
    (
        IOobject::groupName("phaseBlastThermo", phaseName),
        mesh
    ),
    masterName_(masterName),
    phaseName_(phaseName),
    rho_
    (
        mesh.lookupObjectRef<volScalarField>
        (
            IOobject::groupName("rho", phaseName)
        )
    ),
    T_(mesh.lookupObject<volScalarField>(IOobject::groupName("T", masterName_))),
    e_(mesh.lookupObject<volScalarField>(IOobject::groupName("e", masterName_))),
    residualAlpha_("residualAlpha", dimless, 0.0),
    residualRho_("residualRho", dimDensity, 0.0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phaseBlastThermo::~phaseBlastThermo()
{}

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::phaseBlastThermo> Foam::phaseBlastThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName,
    const word& masterName
)
{
    dictionaryConstructorTable::iterator cstrIter =
        blastThermo::lookupCstrIter
        <
            phaseBlastThermo,
            dictionaryConstructorTable
        >
        (
            dict,
            dictionaryConstructorTablePtr_
        );

    return cstrIter()
    (
        mesh,
        dict,
        phaseName,
        masterName
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
