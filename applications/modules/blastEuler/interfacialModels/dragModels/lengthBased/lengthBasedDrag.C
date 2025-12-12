/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
2017-05-18 Jeff Heylmun:    Added support of polydisperse phase models
2025-06-09 Jeff Heylmun:    Added cell based returns
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

#include "lengthBasedDrag.H"
#include "phasePair.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(lengthBased, 0);
    addToRunTimeSelectionTable(dragModel, lengthBased, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::lengthBased::lengthBased
(
    const dictionary& dict,
    const phasePair& pair,
    const bool registerObject
)
:
    dispersedDragModel(dict, pair, registerObject),
    C_("C", dimArea/dimTime, dict)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dragModels::lengthBased::~lengthBased()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::dragModels::lengthBased::CdRe() const
{
    NotImplemented;
    return tmp<volScalarField>();
}


Foam::scalar Foam::dragModels::lengthBased::cellCdRe
(
    const label celli
) const
{
    NotImplemented;
    return 0.0;
}


Foam::tmp<Foam::volScalarField> Foam::dragModels::lengthBased::Ki() const
{
    return
        C_
       *pair_.dispersed().rho()
       *4.0
       /sqr(pair_.dispersed().d());
}


Foam::scalar Foam::dragModels::lengthBased::cellKi
(
    const label celli
) const
{
    return
        C_.value()
       *pair_.dispersed().rho()[celli]
       *4.0
       /sqr(pair_.dispersed().celld(celli));
}

// ************************************************************************* //
