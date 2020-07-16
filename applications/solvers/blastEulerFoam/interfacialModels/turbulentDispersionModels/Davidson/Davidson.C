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

#include "Davidson.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

#include "dragModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulentDispersionModels
{
    defineTypeNameAndDebug(Davidson, 0);
    addToRunTimeSelectionTable
    (
        turbulentDispersionModel,
        Davidson,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Davidson::Davidson
(
    const dictionary& dict,
    const phasePair& pair
)
:
    turbulentDispersionModel(dict, pair),
    residualRe_("residualRe", dimless, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::turbulentDispersionModels::Davidson::~Davidson()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::turbulentDispersionModels::Davidson::D
(
    const label nodei,
    const label nodej
) const
{
    const fvMesh& mesh(pair_.phase1().mesh());
    const volScalarField& alpha1 = pair_.dispersed().volumeFraction(nodei);
    const volScalarField& alpha2 = pair_.continuous().volumeFraction(nodej);
    const dragModel&
        drag
        (
            mesh.lookupObject<dragModel>
            (
                IOobject::groupName(dragModel::typeName, pair_.name())
            )
        );

    volScalarField Cdis
    (
        (4.0/3.0)
       /(
           sqrt(alpha1/max(alpha2, pair_.continuous().residualAlpha()))
          *pair_.continuous().rho()/pair_.dispersed().rho()
         + sqrt(alpha2/max(alpha1, pair_.dispersed().residualAlpha()))
        )
       /(drag.CdRe(nodei, nodej)/max(pair_.Re(nodei, nodej), residualRe_))
    );

    return
        0.75*Cdis
       *pair_.magUr(nodei, nodej)
       *Foam::sqrt(alpha1*alpha2)
       *drag.CdRe(nodei, nodej)
       *pair_.continuous().rho()
       *pair_.continuous().nu()
       /pair_.dispersed().d(nodei)
       /max
        (
            pair_.continuous(),
            pair_.continuous().residualAlpha()
        )*pos0(pair_.dispersed().volumeFraction(nodei) - 0.001);
}


// ************************************************************************* //
