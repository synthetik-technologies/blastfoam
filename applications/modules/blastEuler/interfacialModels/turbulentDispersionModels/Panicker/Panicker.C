/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "Panicker.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

#include "dispersedDragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(Panicker, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        Panicker,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Panicker::Panicker
(
    const dictionary& dict,
    const phasePair& pair
)
:
    dispersedTurbulentDispersionModel(dict, pair),
    Cdis_
    (
        dimensionedScalar::lookupOrDefault
        (
            "Cdis",
            dict,
            dimless,
            4.544
        )
    ),
    residualAlpha_
    (
        "residualAlpha",
        dimless,
        dict.lookupOrDefault<scalar>
        (
            "residualAlpha",
            pair_.dispersed().residualAlpha().value()
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Panicker::~Panicker()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::Panicker::D() const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    tmp<volScalarField> alpha1(pair_.dispersed());
    tmp<volScalarField> d(pair_.dispersed().d());

    const dragModels::dispersedDragModel& drag =
            mesh.lookupObject<dragModels::dispersedDragModel>
            (
                IOobject::groupName(dragModel::typeName, pair_.name())
            );

    scalar b = 0.5;
    scalar a = 1 + b - (1/3);
    return
        0.75
       *drag.CdRe()
       *Cdis_
       *pair_.continuous().rho()
       *sqr(pair_.continuous().nu()/d)
       *pair_.Re()
       *pos0(alpha1() - 0.001)
       *alpha1()*(1 - a*alpha1() + b*sqr(alpha1()));
}


// ************************************************************************* //
